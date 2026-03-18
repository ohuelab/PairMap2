/* ── Pair Mode ─────────────────────────────────────────────────────────────── */

const DEBOUNCE_MS = 350;
let pairCy = null;
let pairInputMode = 'smiles'; // 'smiles' | 'sdf'
let currentSessionId = null;

// BFS-based horizontal layout (source left → target right)
function computeHorizontalLayout(elements) {
  const nodes = elements.filter(e => e.group === 'nodes');
  const edges = elements.filter(e => e.group === 'edges');
  const sourceNode = nodes.find(n => n.data.is_source);
  const sourceId = sourceNode ? sourceNode.data.id : (nodes[0] ? nodes[0].data.id : null);

  // Build adjacency list (undirected)
  const adj = {};
  for (const n of nodes) adj[n.data.id] = [];
  for (const e of edges) {
    adj[e.data.source]?.push(e.data.target);
    adj[e.data.target]?.push(e.data.source);
  }

  // BFS from source
  const dist = {};
  if (sourceId !== null) {
    dist[sourceId] = 0;
    const queue = [sourceId];
    while (queue.length > 0) {
      const cur = queue.shift();
      for (const nb of (adj[cur] || [])) {
        if (!(nb in dist)) {
          dist[nb] = dist[cur] + 1;
          queue.push(nb);
        }
      }
    }
  }
  // Assign distance 0 to any unreached nodes
  for (const n of nodes) {
    if (!(n.data.id in dist)) dist[n.data.id] = 0;
  }

  // Group by rank (distance)
  const layers = {};
  for (const [id, d] of Object.entries(dist)) {
    if (!layers[d]) layers[d] = [];
    layers[d].push(id);
  }

  // Compute positions
  const rankSep = 180;
  const nodeSep = 130;
  const positions = {};
  for (const [rank, ids] of Object.entries(layers)) {
    const x = Number(rank) * rankSep;
    const startY = -((ids.length - 1) * nodeSep) / 2;
    ids.forEach((id, i) => {
      positions[id] = { x, y: startY + i * nodeSep };
    });
  }

  return {
    name: 'preset',
    positions: (node) => positions[node.id()] || { x: 0, y: 0 },
    animate: false,
  };
}

// Debounced SMILES preview
function debounce(fn, ms) {
  let t;
  return (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); };
}

const updatePreviewA = debounce(async (smiles) => {
  await renderMolPreview(smiles, document.getElementById('preview-a'));
}, DEBOUNCE_MS);

const updatePreviewB = debounce(async (smiles) => {
  await renderMolPreview(smiles, document.getElementById('preview-b'));
}, DEBOUNCE_MS);

document.getElementById('smiles-a').addEventListener('input', e => updatePreviewA(e.target.value));
document.getElementById('smiles-b').addEventListener('input', e => updatePreviewB(e.target.value));

/* ── Mode toggle ───────────────────────────────────────────────────────────── */
function setPairMode(mode) {
  pairInputMode = mode;
  const isSMILES = mode === 'smiles';
  document.getElementById('mode-smiles').classList.toggle('active', isSMILES);
  document.getElementById('mode-sdf').classList.toggle('active', !isSMILES);
  document.querySelectorAll('.smiles-mode-inputs').forEach(el => {
    el.style.display = isSMILES ? '' : 'none';
  });
  document.querySelectorAll('.sdf-mode-inputs').forEach(el => {
    el.style.display = isSMILES ? 'none' : '';
  });
  const placeholder = '<div class="mol-preview-placeholder">Enter SMILES to preview</div>';
  document.getElementById('preview-a').innerHTML = placeholder;
  document.getElementById('preview-b').innerHTML = placeholder;
}

document.getElementById('mode-smiles').addEventListener('click', () => setPairMode('smiles'));
document.getElementById('mode-sdf').addEventListener('click', () => setPairMode('sdf'));

/* ── Example button ────────────────────────────────────────────────────────── */
const EXAMPLE_A = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O';  // Ibuprofen
const EXAMPLE_B = 'COc1ccc2cc(ccc2c1)C(C)C(=O)O'; // Naproxen

document.getElementById('pair-example-btn').addEventListener('click', () => {
  setPairMode('smiles');
  document.getElementById('smiles-a').value = EXAMPLE_A;
  document.getElementById('smiles-b').value = EXAMPLE_B;
  updatePreviewA(EXAMPLE_A);
  updatePreviewB(EXAMPLE_B);
});

/* ── SDF file handling ─────────────────────────────────────────────────────── */

/**
 * Render a 3D preview of an SDF file using 3Dmol.js.
 * Falls back to RDKit 2D SVG if 3Dmol is unavailable.
 */
async function render3DPreview(sdfText, container) {
  if (!sdfText) {
    container.innerHTML = '<div class="mol-preview-placeholder">File loaded</div>';
    return;
  }
  if (typeof $3Dmol !== 'undefined') {
    try {
      container.innerHTML = '';
      const viewerDiv = document.createElement('div');
      viewerDiv.style.width = '100%';
      viewerDiv.style.height = '200px';
      viewerDiv.style.position = 'relative';
      container.appendChild(viewerDiv);
      const viewer = $3Dmol.createViewer(viewerDiv, { backgroundColor: 'white' });
      viewer.addModel(sdfText, 'sdf');
      viewer.setStyle({}, { stick: {}, sphere: { scale: 0.25 } });
      viewer.zoomTo();
      viewer.render();
      return;
    } catch (e) {
      console.warn('3Dmol preview failed, falling back to 2D:', e);
    }
  }
  // Fallback: RDKit 2D preview
  await rdkitReady;
  if (!RDKit) {
    container.innerHTML = '<div class="mol-preview-placeholder">File loaded (no preview)</div>';
    return;
  }
  try {
    const end = sdfText.indexOf('$$$$');
    const molblock = end === -1 ? sdfText.trim() : sdfText.slice(0, end).trim();
    const mol = RDKit.get_mol(molblock);
    if (mol && mol.is_valid()) {
      const svg = mol.get_svg(240, 160);
      mol.delete();
      container.innerHTML = svg;
      return;
    }
    if (mol) mol.delete();
  } catch { /* fall through */ }
  container.innerHTML = '<div class="mol-preview-placeholder">File loaded (no preview)</div>';
}

function setupSdfDrop(fileInputId, dropzoneId, filenameId, previewId) {
  const fileInput = document.getElementById(fileInputId);
  const filenameEl = document.getElementById(filenameId);

  async function handleFile(file) {
    if (!file) return;
    filenameEl.textContent = file.name;
    const text = await file.text();
    await render3DPreview(text, document.getElementById(previewId));
  }

  fileInput.addEventListener('change', e => handleFile(e.target.files[0]));

  const dropzone = document.getElementById(dropzoneId);
  dropzone.addEventListener('drop', e => {
    const file = e.dataTransfer?.files[0];
    if (file) {
      const dt = new DataTransfer();
      dt.items.add(file);
      fileInput.files = dt.files;
      handleFile(file);
    }
  });
}

setupSdfDrop('sdf-file-a', 'dropzone-a', 'sdf-filename-a', 'preview-a');
setupSdfDrop('sdf-file-b', 'dropzone-b', 'sdf-filename-b', 'preview-b');

/* ── Collect search/mapgen params ──────────────────────────────────────────── */
function getSearchMapgenParams() {
  const chk = (id) => document.getElementById(id).checked;
  const num = (id) => Number(document.getElementById(id).value);
  const flt = (id) => parseFloat(document.getElementById(id).value);
  return {
    search: {
      is_atom_modfication_enabled: chk('s-atom-mod'),
      cap_ring_with_carbon: chk('s-cap-carbon'),
      cap_ring_with_hydrogen: chk('s-cap-hydrogen'),
      no_backward_search: chk('s-no-backward'),
      use_seed: chk('s-use-seed'),
      ionize: chk('s-ionize'),
      max_intermediate: num('s-max-intermediate'),
    },
    mapgen: {
      maxOptimalPathLength: num('m-max-path'),
      roughScoreThreshold: flt('m-rough-score'),
      minScoreThreshold: flt('m-min-score'),
      CycleLinkThreshold: flt('m-cycle-thresh'),
      optimal_path_mode: chk('m-optimal-mode'),
      squared_sum: chk('m-squared-sum'),
    },
  };
}

/* ── Shared: apply result data to graph + sidebar state ────────────────────── */
async function applyPairResult(data) {
  currentSessionId = data.session_id;

  const graphSection = document.getElementById('pair-graph-section');
  graphSection.classList.add('visible');

  const nodes = data.elements.filter(e => e.group === 'nodes').length;
  const edges = data.elements.filter(e => e.group === 'edges').length;
  document.getElementById('pair-result-stats').textContent =
    `${nodes} nodes · ${edges} edges · ${data.n_intermediates ?? nodes - 2} intermediates`;
  document.getElementById('pair-graph-stats').textContent =
    `${nodes} nodes · ${edges} edges`;

  const layout = computeHorizontalLayout(data.elements);
  pairCy = await renderCytoscape('cy-pair', data, { layout });

  if (pairCy) {
    // Debounced sidebar update based on cytoscape selection state
    let _sidebarTimer = null;
    const updateSidebar = () => {
      clearTimeout(_sidebarTimer);
      _sidebarTimer = setTimeout(() => {
        const selected = pairCy.nodes(':selected');
        if (selected.length === 0) {
          showPairSidebarPlaceholder();
        } else if (selected.length === 1) {
          renderNodeSidebar('pair-sidebar', selected[0].data());
        } else {
          renderOverlaySidebar('pair-sidebar', selected.map(n => n.data()));
        }
      }, 30);
    };

    pairCy.on('tap', 'node', updateSidebar);
    pairCy.on('tap', 'edge', (evt) => {
      const { source, target, similarity } = evt.target.data();
      renderEdgeSidebar('pair-sidebar', source, target, similarity, currentSessionId);
    });
    pairCy.on('tap', (evt) => {
      if (evt.target === pairCy) showPairSidebarPlaceholder();
    });
  }

  // Show Regenerate Map button and download toolbar
  document.getElementById('pair-remap-btn').style.display = '';
  const dlToolbar = document.getElementById('pair-download-toolbar');
  dlToolbar.style.display = 'flex';

  // Wire download buttons
  setupDownloadButtons(currentSessionId);
}

function showPairSidebarPlaceholder() {
  document.getElementById('pair-sidebar').innerHTML = `
    <div class="sidebar-placeholder">
      <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="12" cy="12" r="10"/><line x1="12" y1="8" x2="12" y2="12"/><line x1="12" y1="16" x2="12.01" y2="16"/></svg>
      Click a node to inspect
    </div>`;
}

/* ── Download buttons ──────────────────────────────────────────────────────── */
function setupDownloadButtons(sessionId) {
  const base = API_BASE + `/api/pair/${sessionId}/download`;
  document.getElementById('pair-download-sdf').onclick = () => {
    window.location.href = `${base}/intermediates.sdf`;
  };
  document.getElementById('pair-download-csv').onclick = () => {
    window.location.href = `${base}/links.csv`;
  };
  document.getElementById('pair-download-all-sdf').onclick = () => {
    window.location.href = `${base}/all_intermediates.sdf`;
  };
}

/* ── Generate ──────────────────────────────────────────────────────────────── */
document.getElementById('pair-generate-btn').addEventListener('click', async () => {
  const btn = document.getElementById('pair-generate-btn');
  const spinner = document.getElementById('pair-spinner');
  const spinnerMsg = document.getElementById('pair-spinner-msg');
  const alertEl = document.getElementById('pair-alert');
  const graphSection = document.getElementById('pair-graph-section');

  // Validate before starting
  if (pairInputMode === 'smiles') {
    if (!document.getElementById('smiles-a').value.trim()) {
      showAlert(alertEl, 'Please enter a SMILES for Source.'); return;
    }
    if (!document.getElementById('smiles-b').value.trim()) {
      showAlert(alertEl, 'Please enter a SMILES for Target.'); return;
    }
  } else {
    if (!document.getElementById('sdf-file-a').files[0]) {
      showAlert(alertEl, 'Please provide an SDF file for Source.'); return;
    }
    if (!document.getElementById('sdf-file-b').files[0]) {
      showAlert(alertEl, 'Please provide an SDF file for Target.'); return;
    }
  }

  hideAlert(alertEl);
  btn.disabled = true;
  spinnerMsg.textContent = 'Running PairMap…';
  spinner.classList.add('visible');
  graphSection.classList.remove('visible');

  try {
    let data;
    const { search, mapgen } = getSearchMapgenParams();

    const engine = document.getElementById('pair-engine').value;

    if (pairInputMode === 'smiles') {
      const smilesA = document.getElementById('smiles-a').value.trim();
      const smilesB = document.getElementById('smiles-b').value.trim();
      const res = await fetch(API_BASE + '/api/pair', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles_a: smilesA, smiles_b: smilesB, engine, search, mapgen }),
      });
      if (!res.ok) {
        const err = await res.json().catch(() => ({ detail: res.statusText }));
        throw new Error(typeof err.detail === 'string' ? err.detail : JSON.stringify(err.detail));
      }
      data = await res.json();
    } else {
      const fileA = document.getElementById('sdf-file-a').files[0];
      const fileB = document.getElementById('sdf-file-b').files[0];
      const form = new FormData();
      form.append('file_a', fileA);
      form.append('file_b', fileB);
      form.append('search', JSON.stringify(search));
      form.append('mapgen', JSON.stringify(mapgen));
      form.append('engine', engine);
      const res = await fetch(API_BASE + '/api/pair/sdf', { method: 'POST', body: form });
      if (!res.ok) {
        const err = await res.json().catch(() => ({ detail: res.statusText }));
        throw new Error(typeof err.detail === 'string' ? err.detail : JSON.stringify(err.detail));
      }
      data = await res.json();
    }

    await applyPairResult(data);

  } catch (err) {
    showAlert(alertEl, `Error: ${err.message}`);
  } finally {
    btn.disabled = false;
    spinner.classList.remove('visible');
  }
});

/* ── Regenerate Map ────────────────────────────────────────────────────────── */
document.getElementById('pair-remap-btn').addEventListener('click', async () => {
  if (!currentSessionId) return;

  const btn = document.getElementById('pair-generate-btn');
  const remapBtn = document.getElementById('pair-remap-btn');
  const spinner = document.getElementById('pair-spinner');
  const spinnerMsg = document.getElementById('pair-spinner-msg');
  const alertEl = document.getElementById('pair-alert');

  hideAlert(alertEl);
  btn.disabled = true;
  remapBtn.disabled = true;
  spinnerMsg.textContent = 'Regenerating map…';
  spinner.classList.add('visible');

  try {
    const { mapgen } = getSearchMapgenParams();
    const res = await fetch(API_BASE + '/api/pair/remap', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ session_id: currentSessionId, mapgen }),
    });
    if (!res.ok) {
      const err = await res.json().catch(() => ({ detail: res.statusText }));
      throw new Error(typeof err.detail === 'string' ? err.detail : JSON.stringify(err.detail));
    }
    const data = await res.json();
    await applyPairResult(data);
  } catch (err) {
    showAlert(alertEl, `Error: ${err.message}`);
  } finally {
    btn.disabled = false;
    remapBtn.disabled = false;
    spinner.classList.remove('visible');
  }
});

/* ── Fit button ────────────────────────────────────────────────────────────── */
document.getElementById('pair-fit-btn').addEventListener('click', () => {
  if (pairCy) pairCy.fit(undefined, 40);
});
