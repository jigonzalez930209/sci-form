#!/usr/bin/env node

const fs = require('node:fs');
const path = require('node:path');

const rootDir = path.resolve(__dirname, '..', '..');
const defaultPolicyPath = path.join(rootDir, 'tests', 'fixtures', 'webgpu_regression_policy.json');

function parseArgs(argv) {
  const args = {
    baseline: null,
    candidate: null,
    policy: defaultPolicyPath,
  };

  for (let index = 0; index < argv.length; index += 1) {
    const token = argv[index];
    const next = argv[index + 1];
    switch (token) {
      case '--baseline':
        args.baseline = path.resolve(next);
        index += 1;
        break;
      case '--candidate':
        args.candidate = path.resolve(next);
        index += 1;
        break;
      case '--policy':
        args.policy = path.resolve(next);
        index += 1;
        break;
      default:
        if (token.startsWith('--')) {
          throw new Error(`Unknown option: ${token}`);
        }
    }
  }

  if (!args.baseline || !args.candidate) {
    throw new Error('Usage: node tests/integration/compare_webgpu_reports.js --baseline <file> --candidate <file> [--policy <file>]');
  }

  return args;
}

function readJson(filePath) {
  return JSON.parse(fs.readFileSync(filePath, 'utf8'));
}

function flattenReport(report) {
  const entries = new Map();
  for (const testCase of report.cases) {
    for (const operation of testCase.operations) {
      for (const comparison of operation.comparisons) {
        const key = `${testCase.id}::${operation.id}::${comparison.mode}`;
        entries.set(key, {
          status: comparison.status,
          backend: comparison.backend,
          cpuMeanMs: comparison.cpuMeanMs,
          modeMeanMs: comparison.modeMeanMs,
          speedupVsCpu: comparison.speedupVsCpu,
          maxAbsDiff: comparison.maxAbsDiff,
          rmsDiff: comparison.rmsDiff,
          note: comparison.note,
        });
      }
    }
  }
  return entries;
}

function finite(value) {
  return typeof value === 'number' && Number.isFinite(value);
}

function compareReports(baseline, candidate, policy) {
  const regressions = [];
  const additions = [];
  const baselineEntries = flattenReport(baseline);
  const candidateEntries = flattenReport(candidate);

  for (const [key, baselineEntry] of baselineEntries.entries()) {
    const candidateEntry = candidateEntries.get(key);
    if (!candidateEntry) {
      if (policy.status.failOnMissingComparison) {
        regressions.push(`${key}: missing in candidate report`);
      }
      continue;
    }

    if (
      policy.status.requirePassIfBaselinePassed
      && baselineEntry.status === 'pass'
      && candidateEntry.status !== 'pass'
    ) {
      regressions.push(`${key}: baseline passed but candidate status is ${candidateEntry.status}`);
      continue;
    }

    if (candidateEntry.status !== 'pass') {
      continue;
    }

    const maxAbsLimit = Math.max(
      policy.accuracy.maxAbsFloor,
      (finite(baselineEntry.maxAbsDiff) ? baselineEntry.maxAbsDiff : 0) * policy.accuracy.maxAbsMultiplier
    );
    if (candidateEntry.maxAbsDiff > maxAbsLimit) {
      regressions.push(`${key}: maxAbsDiff ${candidateEntry.maxAbsDiff} exceeded allowed ${maxAbsLimit}`);
    }

    const rmsLimit = Math.max(
      policy.accuracy.rmsFloor,
      (finite(baselineEntry.rmsDiff) ? baselineEntry.rmsDiff : 0) * policy.accuracy.rmsMultiplier
    );
    if (candidateEntry.rmsDiff > rmsLimit) {
      regressions.push(`${key}: rmsDiff ${candidateEntry.rmsDiff} exceeded allowed ${rmsLimit}`);
    }

    if (baselineEntry.speedupVsCpu > 1 && finite(candidateEntry.speedupVsCpu)) {
      const retainedSpeedup = baselineEntry.speedupVsCpu * policy.performance.minSpeedupRetention;
      if (candidateEntry.speedupVsCpu < retainedSpeedup) {
        regressions.push(
          `${key}: speedup ${candidateEntry.speedupVsCpu.toFixed(3)}x fell below retained threshold ${retainedSpeedup.toFixed(3)}x`
        );
      }
    }

    if (
      policy.performance.enforceModeTiming
      && finite(baselineEntry.modeMeanMs)
      && finite(candidateEntry.modeMeanMs)
    ) {
      const modeLimit = baselineEntry.modeMeanMs * policy.performance.maxModeSlowdownRatio;
      if (candidateEntry.modeMeanMs > modeLimit) {
        regressions.push(`${key}: modeMeanMs ${candidateEntry.modeMeanMs.toFixed(3)} exceeded allowed ${modeLimit.toFixed(3)}`);
      }
    }
  }

  for (const key of candidateEntries.keys()) {
    if (!baselineEntries.has(key)) {
      additions.push(key);
    }
  }

  return { regressions, additions };
}

function main() {
  const args = parseArgs(process.argv.slice(2));
  const baseline = readJson(args.baseline);
  const candidate = readJson(args.candidate);
  const policy = readJson(args.policy);
  const { regressions, additions } = compareReports(baseline, candidate, policy);

  console.log(`Baseline:  ${args.baseline}`);
  console.log(`Candidate: ${args.candidate}`);
  console.log(`Policy:    ${args.policy}`);

  if (additions.length > 0) {
    console.log(`New comparisons in candidate: ${additions.join(', ')}`);
  }

  if (regressions.length > 0) {
    console.error('WebGPU regression comparison failed:');
    for (const regression of regressions) {
      console.error(`- ${regression}`);
    }
    process.exit(1);
  }

  console.log('WebGPU regression comparison passed.');
}

main();