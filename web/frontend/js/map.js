/* ── Map Mode ──────────────────────────────────────────────────────────────── */

let mapCy = null;
let mapPollTimer = null;
let mapMolecules = []; // { name, molblock, smiles }

/* ── SDF dropzone / file input ─────────────────────────────────────────────── */
const mapFileInput = document.getElementById('map-sdf-file');
const mapFilenameEl = document.getElementById('map-filename');

mapFileInput.addEventListener('change', async () => {
  const f = mapFileInput.files[0];
  mapFilenameEl.textContent = f ? f.name : '';
  if (f) {
    const text = await f.text();
    await parseSDF(text);
    renderMapMolList();
  }
});

document.getElementById('map-dropzone').addEventListener('drop', async (e) => {
  const file = e.dataTransfer?.files[0];
  if (file) {
    const dt = new DataTransfer();
    dt.items.add(file);
    mapFileInput.files = dt.files;
    mapFilenameEl.textContent = file.name;
    const text = await file.text();
    await parseSDF(text);
    renderMapMolList();
  }
});

document.getElementById('map-example-btn').addEventListener('click', async () => {
  try {
    const res = await fetch(API_BASE + '/examples/p38_example.sdf');
    if (!res.ok) throw new Error('Failed to load example');
    const blob = await res.blob();
    const file = new File([blob], 'p38_example.sdf', { type: 'chemical/x-mdl-sdfile' });
    const dt = new DataTransfer();
    dt.items.add(file);
    mapFileInput.files = dt.files;
    mapFilenameEl.textContent = file.name;
    const text = await blob.text();
    await parseSDF(text);
    renderMapMolList();
  } catch (err) {
    showAlert(document.getElementById('map-alert'), `Error loading example: ${err.message}`);
  }
});

/* ── SDF parser ────────────────────────────────────────────────────────────── */
async function parseSDF(text) {
  mapMolecules = [];
  await rdkitReady;

  const records = text.split('$$$$');
  for (const record of records) {
    const trimmed = record.trim();
    if (!trimmed) continue;

    const molblock = trimmed + '\n$$$$\n';
    const lines = trimmed.split('\n');
    const name = (lines[0] || '').trim() || `mol_${mapMolecules.length + 1}`;

    try {
      const mol = RDKit.get_mol(trimmed);
      if (!mol || !mol.is_valid()) { mol && mol.delete(); continue; }
      const smiles = mol.get_smiles();
      mol.delete();
      mapMolecules.push({ name, molblock, smiles });
    } catch {
      // skip invalid records
    }
  }
}

/* ── Render molecule list ──────────────────────────────────────────────────── */
async function renderMapMolList() {
  const container = document.getElementById('map-mol-list');
  const countEl = document.getElementById('map-mol-count');
  countEl.textContent = mapMolecules.length;

  if (mapMolecules.length === 0) {
    container.innerHTML = '<div class="map-mol-empty">No molecules loaded</div>';
    return;
  }

  await rdkitReady;

  let rows = '';
  for (let i = 0; i < mapMolecules.length; i++) {
    const m = mapMolecules[i];
    let svgHtml = '';
    try {
      const mol = RDKit.get_mol(m.smiles);
      if (mol && mol.is_valid()) {
        const svg = mol.get_svg(80, 60);
        svgHtml = svg;
        mol.delete();
      }
    } catch { /* no thumbnail */ }

    rows += `<tr>
      <td><div class="mol-thumb">${svgHtml}</div></td>
      <td>${escapeHtml(m.name)}</td>
      <td class="mol-smiles-cell" title="${escapeHtml(m.smiles)}">${escapeHtml(m.smiles)}</td>
      <td><button class="btn btn-danger btn-sm" onclick="deleteMapMol(${i})">✕</button></td>
    </tr>`;
  }

  container.innerHTML = `<table>
    <thead><tr><th>Structure</th><th>Name</th><th>SMILES</th><th></th></tr></thead>
    <tbody>${rows}</tbody>
  </table>`;
}

function deleteMapMol(i) {
  mapMolecules.splice(i, 1);
  renderMapMolList();
}

function escapeHtml(str) {
  return str.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}

/* ── Add SMILES manually ───────────────────────────────────────────────────── */
document.getElementById('map-add-smiles-btn').addEventListener('click', async () => {
  const input = document.getElementById('map-smiles-input');
  const smiles = input.value.trim();
  if (!smiles) return;

  const alertEl = document.getElementById('map-alert');
  await rdkitReady;

  try {
    const mol = RDKit.get_mol(smiles);
    if (!mol || !mol.is_valid()) { mol && mol.delete(); showAlert(alertEl, 'Invalid SMILES'); return; }
    const molblock = mol.get_molblock() + '\n$$$$\n';
    const name = `mol_${mapMolecules.length + 1}`;
    mol.delete();
    mapMolecules.push({ name, molblock, smiles });
    input.value = '';
    hideAlert(alertEl);
    renderMapMolList();
  } catch (err) {
    showAlert(alertEl, `Invalid SMILES: ${err.message}`);
  }
});

/* ── Clear All ─────────────────────────────────────────────────────────────── */
document.getElementById('map-clear-all-btn').addEventListener('click', () => {
  mapMolecules = [];
  mapFileInput.value = '';
  mapFilenameEl.textContent = '';
  renderMapMolList();
});

/* ── Collect map config ────────────────────────────────────────────────────── */
function getMapConfig() {
  const num = (id) => Number(document.getElementById(id).value);
  const flt = (id) => parseFloat(document.getElementById(id).value);
  const chk = (id) => document.getElementById(id).checked;

  return {
    similarity_threshold: flt('map-sim-thresh'),
    max_intermediate: num('map-max-inter'),
    max: num('map-max-degree'),
    max_dist_from_actives: num('map-max-dist'),
    max_path_length: num('map-max-path'),
    allow_tree: chk('map-allow-tree'),
    radial: chk('map-radial'),
    ionize: chk('map-ionize'),
    maxOptimalPathLength: num('map-max-opt-path'),
    roughScoreThreshold: flt('map-rough-score'),
    minScoreThreshold: flt('map-min-score'),
    optimal_path_mode: chk('map-optimal-mode'),
  };
}

/* ── Submit job ────────────────────────────────────────────────────────────── */
document.getElementById('map-submit-btn').addEventListener('click', async () => {
  const alertEl = document.getElementById('map-alert');

  if (mapMolecules.length === 0) {
    showAlert(alertEl, 'Please load an SDF file or add molecules.');
    return;
  }

  hideAlert(alertEl);

  const engine = document.getElementById('map-engine').value;
  const config = getMapConfig();

  // Reconstruct SDF from mapMolecules
  const sdfText = mapMolecules.map(m => m.molblock).join('\n');
  const sdfBlob = new Blob([sdfText], { type: 'chemical/x-mdl-sdfile' });
  const sdfFile = new File([sdfBlob], 'input.sdf', { type: 'chemical/x-mdl-sdfile' });

  const fd = new FormData();
  fd.append('file', sdfFile);
  fd.append('engine', engine);
  fd.append('config', JSON.stringify(config));

  try {
    const res = await apiFetch(API_BASE + '/api/map/jobs', { method: 'POST', body: fd });
    if (!res.ok) await throwApiError(res);
    const job = await res.json();
    showAlert(alertEl, `Job submitted: ${job.id}`, 'info');
    refreshMapJobs();
    startMapPolling(job.id);
  } catch (err) {
    showAlert(alertEl, `Error: ${err.message}`);
  }
});

/* ── Job list ──────────────────────────────────────────────────────────────── */
document.getElementById('map-refresh-btn').addEventListener('click', refreshMapJobs);

async function refreshMapJobs() {
  const tbody = document.getElementById('map-jobs-body');
  try {
    const res = await apiFetch(API_BASE + '/api/map/jobs');
    const data = await res.json();
    const jobs = data.jobs || [];

    if (jobs.length === 0) {
      tbody.innerHTML = '<tr><td colspan="6" style="text-align:center;color:var(--text-muted);padding:2rem;">No jobs yet</td></tr>';
      return;
    }

    tbody.innerHTML = jobs.map(job => {
      const created = new Date(job.created_at).toLocaleString();
      const statusClass = `status-${job.status}`;
      const actions = [
        job.status === 'completed'
          ? `<button class="btn btn-secondary btn-sm" onclick="viewMapGraph('${job.id}')">View</button>
             <a class="btn btn-secondary btn-sm" href="${API_BASE}/api/map/jobs/${job.id}/artifacts/intermediate_mols.sdf" download>SDF</a>`
          : '',
        ['queued', 'running'].includes(job.status)
          ? `<button class="btn btn-danger btn-sm" onclick="cancelMapJob('${job.id}')">Cancel</button>`
          : '',
      ].join(' ');

      return `<tr>
        <td class="job-id-cell" title="${job.id}">${job.id.slice(0, 8)}…</td>
        <td>${job.engine}</td>
        <td><span class="status-badge ${statusClass}">${job.status}</span></td>
        <td style="color:var(--text-secondary);font-size:12px;">${created}</td>
        <td style="font-size:12px;color:var(--text-muted);">${job.progress || '—'}</td>
        <td style="display:flex;gap:0.4rem;flex-wrap:wrap;">${actions}</td>
      </tr>`;
    }).join('');

  } catch (err) {
    console.error('Failed to refresh map jobs:', err);
  }
}

/* ── Polling ───────────────────────────────────────────────────────────────── */
function startMapPolling(jobId) {
  if (mapPollTimer) clearInterval(mapPollTimer);
  mapPollTimer = setInterval(async () => {
    try {
      const res = await apiFetch(API_BASE + `/api/map/jobs/${jobId}`);
      const job = await res.json();
      refreshMapJobs();
      if (!['queued', 'running'].includes(job.status)) {
        clearInterval(mapPollTimer);
        mapPollTimer = null;
      }
    } catch {
      clearInterval(mapPollTimer);
      mapPollTimer = null;
    }
  }, 3000);
}

/* ── Cancel job ────────────────────────────────────────────────────────────── */
async function cancelMapJob(jobId) {
  try {
    await apiFetch(API_BASE + `/api/map/jobs/${jobId}/cancel`, { method: 'POST' });
    refreshMapJobs();
  } catch (err) {
    console.error('Cancel failed:', err);
  }
}

/* ── View graph ────────────────────────────────────────────────────────────── */
async function viewMapGraph(jobId) {
  const graphSection = document.getElementById('map-graph-section');
  document.getElementById('map-graph-job-label').textContent = jobId.slice(0, 8) + '…';

  try {
    const res = await apiFetch(API_BASE + `/api/map/jobs/${jobId}/graph`);
    if (!res.ok) { showAlert(document.getElementById('map-alert'), 'Graph not available yet.'); return; }
    const data = await res.json();

    graphSection.classList.add('visible');
    graphSection.scrollIntoView({ behavior: 'smooth', block: 'start' });

    const toolbar = document.getElementById('map-download-toolbar');
    toolbar.style.display = 'flex';

    document.getElementById('map-download-sdf').onclick = () => {
      window.location.href = API_BASE + `/api/map/jobs/${jobId}/artifacts/intermediate_mols.sdf`;
    };

    document.getElementById('map-download-links-csv').onclick = () => {
      const edges = data.elements.filter(e => e.group === 'edges');
      const rows = ['source,target,similarity,bad_edge'];
      for (const e of edges) {
        const { source, target, similarity, bad_edge } = e.data;
        rows.push(`${source},${target},${similarity != null ? similarity : ''},${bad_edge ? 'true' : 'false'}`);
      }
      const blob = new Blob([rows.join('\n')], { type: 'text/csv' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = `map_links_${jobId.slice(0, 8)}.csv`;
      a.click(); URL.revokeObjectURL(url);
    };

    document.getElementById('map-download-json').onclick = () => {
      const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = `map_graph_${jobId.slice(0, 8)}.json`;
      a.click(); URL.revokeObjectURL(url);
    };

    const nodes = data.elements.filter(e => e.group === 'nodes').length;
    const edges = data.elements.filter(e => e.group === 'edges').length;
    document.getElementById('map-graph-stats').textContent = `${nodes} nodes · ${edges} edges`;

    mapCy = await renderCytoscape('cy-map', data);

    if (mapCy) {
      mapCy.on('tap', 'node', async (evt) => {
        await renderNodeSidebar('map-sidebar', evt.target.data());
      });
      mapCy.on('tap', 'edge', async (evt) => {
        await renderMapEdgeSidebar(jobId, evt.target.data());
      });
      mapCy.on('tap', (evt) => {
        if (evt.target === mapCy) {
          document.getElementById('map-sidebar').innerHTML = `
            <div class="sidebar-placeholder">
              <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="12" cy="12" r="10"/><line x1="12" y1="8" x2="12" y2="12"/><line x1="12" y1="16" x2="12.01" y2="16"/></svg>
              Click a node or edge to inspect
            </div>`;
        }
      });
    }

  } catch (err) {
    showAlert(document.getElementById('map-alert'), `Error loading graph: ${err.message}`);
  }
}

/* ── Map edge MCS sidebar ──────────────────────────────────────────────────── */
async function renderMapEdgeSidebar(jobId, edgeData) {
  const sidebar = document.getElementById('map-sidebar');
  const { source, target, similarity, bad_edge } = edgeData;

  sidebar.innerHTML = `
    <div class="sidebar-placeholder">
      <div class="spinner" style="margin:auto"></div>
      <span style="font-size:12px;color:var(--text-muted);margin-top:0.5rem;">Computing MCS…</span>
    </div>`;

  try {
    const res = await apiFetch(API_BASE + `/api/map/jobs/${jobId}/mcs/${source}/${target}`);
    if (!res.ok) throw new Error(await res.text());
    const data = await res.json();

    const simText = similarity != null ? Number(similarity).toFixed(3) : '—';
    const badBadge = bad_edge ? `<span class="status-badge status-failed" style="font-size:10px;">bad edge</span>` : '';

    sidebar.innerHTML = `
      <div class="edge-mcs-info">
        <div class="node-name" style="font-size:13px;">Edge: ${data.label_a} → ${data.label_b}</div>
        <div style="font-size:12px;color:var(--text-muted);display:flex;align-items:center;gap:0.5rem;">
          Similarity: ${simText} ${badBadge}
        </div>
        <div class="mcs-legend">
          <span class="mcs-legend-item mcs-common">■ Common (${data.n_common})</span>
          <span class="mcs-legend-item mcs-deleted">■ Deleted (${data.n_deleted})</span>
          <span class="mcs-legend-item mcs-inserted">■ Inserted (${data.n_inserted})</span>
        </div>
        <div class="mcs-mol-pair">
          <div class="mcs-mol-item">
            <div class="mcs-mol-label">${data.label_a}</div>
            <div class="mcs-mol-svg">${data.svg_a}</div>
          </div>
          <div class="mcs-mol-item">
            <div class="mcs-mol-label">${data.label_b}</div>
            <div class="mcs-mol-svg">${data.svg_b}</div>
          </div>
        </div>
      </div>`;
  } catch (e) {
    const placeholder = document.createElement('div');
    placeholder.className = 'sidebar-placeholder';
    placeholder.style.color = 'var(--red)';
    placeholder.textContent = `MCS computation failed: ${e.message}`;
    sidebar.innerHTML = '';
    sidebar.appendChild(placeholder);
  }
}

/* ── Fit button ────────────────────────────────────────────────────────────── */
document.getElementById('map-fit-btn').addEventListener('click', () => {
  if (mapCy) mapCy.fit(undefined, 40);
});

/* ── Reset layout button ───────────────────────────────────────────────────── */
document.getElementById('map-reset-btn').addEventListener('click', () => {
  if (mapCy) resetPositions(mapCy);
});

/* ── Init ──────────────────────────────────────────────────────────────────── */
renderMapMolList();
refreshMapJobs();
