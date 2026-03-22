import init, {
  embed,
  initThreadPool,
  init_webgpu,
  compute_esp_grid_accelerated_typed,
  compute_esp_grid_info,
  eht_orbital_grid_accelerated_typed,
  eht_density_grid_accelerated_typed,
  eht_volumetric_grid_info,
  compute_orbital_mesh_accelerated,
} from '../../../pkg/sci_form_wasm.js';

const output = document.getElementById('output');

function print(value) {
  output.textContent = `${output.textContent}\n${value}`.trim();
}

function webgpuHelpMessage(error) {
  const lines = [
    'WebGPU is required for this example.',
    'Enable WebGPU in a supported browser and reload the page.',
    'Suggested browsers/settings:',
    '- Chrome / Edge: open chrome://flags/#enable-unsafe-webgpu and enable WebGPU',
    '- Firefox Nightly: enable dom.webgpu.enabled in about:config',
    '- Make sure hardware acceleration is enabled in the browser settings',
  ];

  if (error) {
    lines.push(`Runtime error: ${error}`);
  }

  return lines.join('\n');
}

async function main() {
  await init();
  await initThreadPool(Math.max(1, navigator.hardwareConcurrency ?? 4));

  if (!navigator.gpu) {
    print(webgpuHelpMessage('navigator.gpu is not available in this browser'));
    return;
  }

  const gpuStatus = JSON.parse(await init_webgpu());
  if (!gpuStatus.initialized) {
    print(webgpuHelpMessage(gpuStatus.error));
    return;
  }

  print(`WebGPU ready: ${gpuStatus.backend} | max storage buffer ${gpuStatus.max_storage_buffer_size}`);

  const ethanol = JSON.parse(embed('CCO', 42));
  const elements = JSON.stringify(ethanol.elements);
  const coords = JSON.stringify(ethanol.coords);

  const espInfo = JSON.parse(compute_esp_grid_info(elements, coords, 0.4, 3.0));
  const espValues = await compute_esp_grid_accelerated_typed(
    elements,
    coords,
    0.4,
    3.0,
    'hybrid'
  );
  print(`ESP typed grid: len=${espValues.length} | dims=${espInfo.dims.join('x')}`);

  const volumetricInfo = JSON.parse(eht_volumetric_grid_info(elements, coords, 0.4, 3.0));
  const orbitalValues = await eht_orbital_grid_accelerated_typed(
    elements,
    coords,
    0,
    0.4,
    3.0,
    'gpu'
  );
  print(`Orbital typed grid: len=${orbitalValues.length} | dims=${volumetricInfo.dims.join('x')}`);

  const densityValues = await eht_density_grid_accelerated_typed(
    elements,
    coords,
    0.4,
    3.0,
    'auto'
  );
  print(`Density typed grid: len=${densityValues.length}`);

  const mesh = JSON.parse(
    await compute_orbital_mesh_accelerated(
      elements,
      coords,
      'eht',
      0,
      0.4,
      3.0,
      0.02,
      'hybrid'
    )
  );
  print(`Mesh backend: ${mesh.backend} | triangles=${mesh.mesh.num_triangles}`);
}

main().catch((error) => {
  output.textContent = `Error:\n${error?.stack ?? error}`;
});