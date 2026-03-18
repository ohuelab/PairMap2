/* ── Map Mode ──────────────────────────────────────────────────────────────── */

let mapCy = null;
let mapPollTimer = null;

/* ── SDF dropzone / file input ─────────────────────────────────────────────── */
const mapFileInput = document.getElementById('map-sdf-file');
const mapFilenameEl = document.getElementById('map-filename');

mapFileInput.addEventListener('change', () => {
  const f = mapFileInput.files[0];
  mapFilenameEl.textContent = f ? f.name : '';
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
    maxOptimalPathLength: num('map-max-opt-path'),
    roughScoreThreshold: flt('map-rough-score'),
    minScoreThreshold: flt('map-min-score'),
    optimal_path_mode: chk('map-optimal-mode'),
  };
}

/* ── Submit job ────────────────────────────────────────────────────────────── */
document.getElementById('map-submit-btn').addEventListener('click', async () => {
  const alertEl = document.getElementById('map-alert');
  const file = mapFileInput.files[0];
  if (!file) { showAlert(alertEl, 'Please select an SDF file.'); return; }

  hideAlert(alertEl);

  const engine = document.getElementById('map-engine').value;
  const config = getMapConfig();

  const fd = new FormData();
  fd.append('file', file);
  fd.append('engine', engine);
  fd.append('config', JSON.stringify(config));

  try {
    const res = await fetch(API_BASE + '/api/map/jobs', { method: 'POST', body: fd });
    if (!res.ok) {
      const err = await res.json().catch(() => ({ detail: res.statusText }));
      throw new Error(typeof err.detail === 'string' ? err.detail : JSON.stringify(err.detail));
    }
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
    const res = await fetch(API_BASE + '/api/map/jobs');
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
      const res = await fetch(API_BASE + `/api/map/jobs/${jobId}`);
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
    await fetch(API_BASE + `/api/map/jobs/${jobId}/cancel`, { method: 'POST' });
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
    const res = await fetch(API_BASE + `/api/map/jobs/${jobId}/graph`);
    if (!res.ok) { showAlert(document.getElementById('map-alert'), 'Graph not available yet.'); return; }
    const data = await res.json();

    graphSection.classList.add('visible');
    graphSection.scrollIntoView({ behavior: 'smooth', block: 'start' });

    const nodes = data.elements.filter(e => e.group === 'nodes').length;
    const edges = data.elements.filter(e => e.group === 'edges').length;
    document.getElementById('map-graph-stats').textContent = `${nodes} nodes · ${edges} edges`;

    mapCy = await renderCytoscape('cy-map', data);

    if (mapCy) {
      mapCy.on('tap', 'node', async (evt) => {
        await renderNodeSidebar('map-sidebar', evt.target.data());
      });
      mapCy.on('tap', (evt) => {
        if (evt.target === mapCy) {
          document.getElementById('map-sidebar').innerHTML = `
            <div class="sidebar-placeholder">
              <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="12" cy="12" r="10"/><line x1="12" y1="8" x2="12" y2="12"/><line x1="12" y1="16" x2="12.01" y2="16"/></svg>
              Click a node to inspect
            </div>`;
        }
      });
    }

  } catch (err) {
    showAlert(document.getElementById('map-alert'), `Error loading graph: ${err.message}`);
  }
}

/* ── Fit button ────────────────────────────────────────────────────────────── */
document.getElementById('map-fit-btn').addEventListener('click', () => {
  if (mapCy) mapCy.fit(undefined, 40);
});

/* ── Init ──────────────────────────────────────────────────────────────────── */
refreshMapJobs();
