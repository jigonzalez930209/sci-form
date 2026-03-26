#!/usr/bin/env node

const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const { spawn } = require('node:child_process');
const { chromium } = require('playwright');

const rootDir = path.resolve(__dirname, '..', '..');
const harnessDir = path.join(rootDir, 'crates', 'wasm', 'examples', 'vite-webgpu-hybrid');
const defaultReportPath = path.join(rootDir, 'build', 'webgpu', 'latest-report.json');

function parseArgs(argv) {
  const args = {
    port: 4173,
    iterations: 3,
    timeoutMs: 180000,
    profile: 'smoke',
    caseId: null,
    reportOut: defaultReportPath,
    allowSkips: false,
    allowNoWebgpu: false,
    url: null,
  };

  for (let index = 0; index < argv.length; index += 1) {
    const token = argv[index];
    const next = argv[index + 1];
    switch (token) {
      case '--port':
        args.port = Number.parseInt(next, 10);
        index += 1;
        break;
      case '--iterations':
        args.iterations = Number.parseInt(next, 10);
        index += 1;
        break;
      case '--timeout-ms':
        args.timeoutMs = Number.parseInt(next, 10);
        index += 1;
        break;
      case '--profile':
        args.profile = next;
        index += 1;
        break;
      case '--case':
        args.caseId = next;
        index += 1;
        break;
      case '--report-out':
        args.reportOut = path.resolve(next);
        index += 1;
        break;
      case '--url':
        args.url = next;
        index += 1;
        break;
      case '--allow-skips':
        args.allowSkips = true;
        break;
      case '--allow-no-webgpu':
        args.allowNoWebgpu = true;
        break;
      default:
        if (token.startsWith('--')) {
          throw new Error(`Unknown option: ${token}`);
        }
    }
  }

  return args;
}

function npmCommand() {
  return process.platform === 'win32' ? 'npm.cmd' : 'npm';
}

function ensureDirectory(filePath) {
  fs.mkdirSync(path.dirname(filePath), { recursive: true });
}

function httpGet(url) {
  return new Promise((resolve, reject) => {
    const req = http.get(url, (res) => {
      res.resume();
      resolve(res.statusCode || 0);
    });
    req.on('error', reject);
  });
}

async function waitForServer(url, timeoutMs) {
  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {
    try {
      const status = await httpGet(url);
      if (status >= 200 && status < 500) {
        return;
      }
    } catch (_error) {
    }
    await new Promise((resolve) => setTimeout(resolve, 500));
  }
  throw new Error(`Timed out waiting for preview server at ${url}`);
}

function startPreviewServer(port) {
  const child = spawn(
    npmCommand(),
    ['--prefix', harnessDir, 'run', 'preview', '--', '--host', '127.0.0.1', '--port', String(port), '--strictPort'],
    {
      cwd: rootDir,
      stdio: ['ignore', 'pipe', 'pipe'],
    }
  );

  child.stdout.on('data', (chunk) => process.stdout.write(chunk));
  child.stderr.on('data', (chunk) => process.stderr.write(chunk));

  return child;
}

function buildHarness() {
  return new Promise((resolve, reject) => {
    const child = spawn(
      npmCommand(),
      ['--prefix', harnessDir, 'run', 'build'],
      {
        cwd: rootDir,
        stdio: ['ignore', 'pipe', 'pipe'],
      }
    );

    child.stdout.on('data', (chunk) => process.stdout.write(chunk));
    child.stderr.on('data', (chunk) => process.stderr.write(chunk));
    child.on('error', reject);
    child.on('exit', (code) => {
      if (code === 0) {
        resolve();
      } else {
        reject(new Error(`Harness build failed with exit code ${code}`));
      }
    });
  });
}

function collectStatuses(report) {
  const failures = [];
  const skips = [];
  for (const testCase of report.cases) {
    for (const operation of testCase.operations) {
      for (const comparison of operation.comparisons) {
        const label = `${testCase.id}/${operation.id}/${comparison.mode}`;
        if (comparison.status === 'fail') {
          failures.push(`${label}: ${comparison.note}`);
        } else if (comparison.status === 'skip') {
          skips.push(`${label}: ${comparison.note}`);
        }
      }
    }
  }
  return { failures, skips };
}

async function collectHarnessReport(url, timeoutMs) {
  const browser = await chromium.launch({
    headless: true,
    args: [
      '--enable-unsafe-webgpu',
      '--ignore-gpu-blocklist',
      '--enable-features=Vulkan,UseSkiaRenderer',
    ],
  });

  const page = await browser.newPage();
  const consoleLines = [];
  page.on('console', (message) => {
    consoleLines.push(`[${message.type()}] ${message.text()}`);
  });

  try {
    await page.goto(url, { waitUntil: 'domcontentloaded', timeout: timeoutMs });
    await page.waitForFunction(
      () => Boolean(window.__SCI_FORM_WEBGPU_REPORT__ && window.__SCI_FORM_WEBGPU_REPORT__.completed),
      { timeout: timeoutMs }
    );
    const report = await page.evaluate(() => window.__SCI_FORM_WEBGPU_REPORT__);
    return { report, consoleLines };
  } finally {
    await page.close();
    await browser.close();
  }
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  const harnessUrl = args.url
    || `http://127.0.0.1:${args.port}/benchmark.html?profile=${encodeURIComponent(args.profile)}&iterations=${args.iterations}${args.caseId ? `&case=${encodeURIComponent(args.caseId)}` : ''}`;

  let preview = null;
  if (!args.url) {
    await buildHarness();
    preview = startPreviewServer(args.port);
    await waitForServer(`http://127.0.0.1:${args.port}/benchmark.html`, args.timeoutMs);
  }

  try {
    const { report, consoleLines } = await collectHarnessReport(harnessUrl, args.timeoutMs);
    report.runner = {
      harnessUrl,
      collectedAt: new Date().toISOString(),
      consoleLines,
    };

    ensureDirectory(args.reportOut);
    fs.writeFileSync(args.reportOut, `${JSON.stringify(report, null, 2)}\n`);

    const { failures, skips } = collectStatuses(report);
    console.log(`Saved report to ${args.reportOut}`);
    console.log(`Comparisons: ${report.summary.passed}/${report.summary.completed} passed, ${report.summary.skipped} skipped`);

    if (!report.webgpu.initialized && !args.allowNoWebgpu) {
      throw new Error(`WebGPU was not initialized: ${report.webgpu.error || 'unknown error'}`);
    }
    if (failures.length > 0) {
      throw new Error(`Harness reported failing comparisons:\n${failures.join('\n')}`);
    }
    if (skips.length > 0 && !args.allowSkips) {
      throw new Error(`Harness reported skipped comparisons:\n${skips.join('\n')}`);
    }
  } finally {
    if (preview) {
      preview.kill('SIGTERM');
    }
  }
}

main().catch((error) => {
  console.error(error?.stack ?? String(error));
  process.exit(1);
});