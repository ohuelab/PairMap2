/* ── API base URL ───────────────────────────────────────────────────────────── */
// In production (Cloudflare Pages), point to the Render.com backend.
// In development (localhost), use the same origin (empty string).
const API_BASE = window.location.hostname === 'localhost'
  ? ''
  : 'https://pairmap-api.onrender.com';

/* ── RDKit.js init ─────────────────────────────────────────────────────────── */
let RDKit = null;
const rdkitReady = (async () => {
  try {
    const mod = await window.initRDKitModule({
      locateFile: (f) =>
        `https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/${f}`,
    });
    RDKit = mod;
    return mod;
  } catch (e) {
    console.warn('RDKit.js failed to load:', e);
    return null;
  }
})();

/**
 * Render a SMILES string to an SVG string.
 * Returns null on failure.
 */
async function smilesTo2DSVG(smiles, w = 240, h = 160) {
  await rdkitReady;
  if (!RDKit || !smiles) return null;
  try {
    const mol = RDKit.get_mol(smiles);
    if (!mol.is_valid()) { mol.delete(); return null; }
    const svg = mol.get_svg(w, h);
    mol.delete();
    return svg;
  } catch {
    return null;
  }
}

/**
 * Convert an SVG string to a data: URL.
 */
function svgToDataUrl(svg) {
  return 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svg)));
}

/**
 * Render mol to a DOM container element.
 * Shows placeholder or error message when SMILES is invalid.
 */
async function renderMolPreview(smiles, container) {
  if (!smiles || !smiles.trim()) {
    container.innerHTML = '<div class="mol-preview-placeholder">Enter SMILES to preview</div>';
    return;
  }
  const svg = await smilesTo2DSVG(smiles.trim());
  if (svg) {
    container.innerHTML = svg;
  } else {
    container.innerHTML = '<div class="mol-preview-error">Invalid SMILES</div>';
  }
}

/* ── Tab switching ─────────────────────────────────────────────────────────── */
document.querySelectorAll('.tab-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    const tab = btn.dataset.tab;
    document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
    document.querySelectorAll('.tab-content').forEach(s => s.classList.remove('active'));
    btn.classList.add('active');
    document.getElementById(`tab-${tab}`).classList.add('active');
  });
});

/* ── Collapsible param sections ────────────────────────────────────────────── */
document.querySelectorAll('.param-section-header').forEach(header => {
  const sectionId = header.dataset.section;
  const body = document.getElementById(sectionId);
  header.addEventListener('click', () => {
    const isOpen = header.classList.toggle('open');
    body.classList.toggle('open', isOpen);
  });
});

/* ── Slider sync: update display value on input ────────────────────────────── */
document.querySelectorAll('input[type="range"]').forEach(slider => {
  const valEl = document.getElementById(slider.id + '-val');
  if (!valEl) return;
  const update = () => { valEl.textContent = parseFloat(slider.value).toFixed(2); };
  slider.addEventListener('input', update);
  update();
});

/* ── Dropzone drag-over highlight ──────────────────────────────────────────── */
document.querySelectorAll('.dropzone').forEach(dz => {
  dz.addEventListener('dragover', e => { e.preventDefault(); dz.classList.add('drag-over'); });
  dz.addEventListener('dragleave', () => dz.classList.remove('drag-over'));
  dz.addEventListener('drop', e => { e.preventDefault(); dz.classList.remove('drag-over'); });
});

/* ── Show / hide alert ─────────────────────────────────────────────────────── */
function showAlert(el, msg, type = 'error') {
  el.textContent = msg;
  el.className = `alert visible ${type}`;
}

function hideAlert(el) {
  el.className = 'alert';
  el.textContent = '';
}

/* ── Cytoscape helpers ─────────────────────────────────────────────────────── */

/**
 * Build Cytoscape.js stylesheet for a light molecular graph.
 */
function makeCyStyle() {
  return [
    {
      selector: 'node',
      style: {
        'width': 72,
        'height': 72,
        'background-color': '#ffffff',
        'background-fit': 'contain',
        'background-clip': 'none',
        'background-image': 'none',
        'border-width': 3,
        'border-color': '#0969da',
        'label': 'data(label)',
        'text-valign': 'bottom',
        'text-halign': 'center',
        'font-size': 10,
        'color': '#656d76',
        'text-margin-y': 10,
        'text-max-width': '90px',
        'text-overflow-wrap': 'ellipsis',
        'text-wrap': 'ellipsis',
      },
    },
    {
      selector: 'node[?is_source]',
      style: { 'border-color': '#0969da', 'border-width': 3 },
    },
    {
      selector: 'node[?is_target]',
      style: { 'border-color': '#1a7f37', 'border-width': 3 },
    },
    {
      selector: 'node[?intermediate]',
      style: { 'border-color': '#9a6700', 'border-width': 2 },
    },
    {
      selector: 'node[?active][!intermediate]',
      style: { 'border-color': '#0969da', 'border-width': 3 },
    },
    {
      selector: 'node:selected',
      style: {
        'border-color': '#e36209',
        'border-width': 4,
        'background-color': '#fff8f0',
      },
    },
    {
      selector: 'edge',
      style: {
        'width': (ele) => Math.max(1.5, (ele.data('similarity') || 0) * 5),
        'line-color': (ele) => ele.data('bad_edge') ? '#cf222e' : '#b0b8c1',
        'label': (ele) => (ele.data('similarity') || 0).toFixed(2),
        'font-size': 9,
        'color': '#9198a1',
        'text-background-color': '#ffffff',
        'text-background-opacity': 1,
        'text-background-padding': '3px',
        'curve-style': 'bezier',
        'opacity': 0.85,
      },
    },
    {
      selector: 'edge:selected',
      style: { 'line-color': '#e36209', 'width': 3, 'opacity': 1 },
    },
  ];
}

/**
 * Create or reinitialise a Cytoscape instance.
 *
 * @param {string} containerId - DOM element id for the cytoscape container.
 * @param {object} data - Response data containing `elements` array.
 * @param {object} [options] - Optional overrides. Supports `options.layout`.
 */
async function renderCytoscape(containerId, data, options = {}) {
  const container = document.getElementById(containerId);
  if (!container) return null;

  // Destroy existing instance
  if (window._cyInstances && window._cyInstances[containerId]) {
    window._cyInstances[containerId].destroy();
  }
  if (!window._cyInstances) window._cyInstances = {};

  const nodeCount = (data.elements || []).filter(e => e.group === 'nodes').length;

  const defaultLayout = {
    name: 'cose',
    padding: 48,
    animate: nodeCount <= 40,
    animationDuration: 400,
    randomize: false,
    nodeRepulsion: () => 10000,
    idealEdgeLength: () => 100,
    edgeElasticity: () => 100,
    nestingFactor: 1.2,
    gravity: 0.25,
    numIter: 1000,
    coolingFactor: 0.99,
    minTemp: 1.0,
  };

  const layoutConfig = options.layout || defaultLayout;

  const cy = cytoscape({
    container,
    elements: data.elements || [],
    style: makeCyStyle(),
    wheelSensitivity: 0.3,
  });

  window._cyInstances[containerId] = cy;

  const layoutInstance = cy.layout(layoutConfig);

  // Load molecule images after layout completes (works for both animated and instant layouts)
  const loadImages = async () => {
    for (const node of cy.nodes()) {
      // Prefer pre-computed MCS-aligned SVG from backend
      const alignedSvg = node.data('aligned_svg');
      if (alignedSvg) {
        node.style({ 'background-image': svgToDataUrl(alignedSvg), 'background-color': '#ffffff' });
        continue;
      }
      // Fallback: render SMILES client-side via RDKit.js
      await rdkitReady;
      if (!RDKit) continue;
      const smiles = node.data('smiles');
      if (!smiles) continue;
      const svg = await smilesTo2DSVG(smiles, 100, 80);
      if (svg) {
        node.style({
          'background-image': svgToDataUrl(svg),
          'background-color': '#ffffff',
        });
      }
    }
  };

  if (layoutConfig.animate) {
    layoutInstance.one('layoutstop', loadImages);
  } else {
    layoutInstance.run();
    loadImages();
    return cy;
  }

  layoutInstance.run();

  return cy;
}

/* ── 3D molecule viewer helpers ────────────────────────────────────────────── */

/**
 * Create a 3Dmol.js viewer in a container div.
 * Adds stick + small sphere style and renders.
 */
function create3DViewer(containerEl, molblock, color) {
  if (typeof $3Dmol === 'undefined') return null;
  try {
    containerEl.style.width = '100%';
    containerEl.style.height = '220px';
    containerEl.style.position = 'relative';
    const viewer = $3Dmol.createViewer(containerEl, { backgroundColor: 'white' });
    viewer.addModel(molblock, 'sdf');
    const styleSpec = { stick: {}, sphere: { scale: 0.25 } };
    if (color) {
      styleSpec.stick.color = color;
      styleSpec.sphere.color = color;
    }
    viewer.setStyle({}, styleSpec);
    viewer.zoomTo();
    viewer.render();
    return viewer;
  } catch (e) {
    console.warn('3Dmol viewer error:', e);
    return null;
  }
}

/**
 * Render the node info sidebar (3D if molblock present, 2D SVG fallback).
 * MW/LogP/HBA/HBD properties are intentionally omitted.
 */
async function renderNodeSidebar(sidebarId, nodeData) {
  const sidebar = document.getElementById(sidebarId);
  if (!sidebar) return;

  const { id, label, smiles, is_source, is_target, molblock } = nodeData;

  let typeClass = 'intermediate';
  let typeLabel = 'Intermediate';
  if (is_source) { typeClass = 'source'; typeLabel = 'Source'; }
  else if (is_target) { typeClass = 'target'; typeLabel = 'Target'; }

  const viewerId = `sidebar-3d-${Date.now()}`;
  sidebar.innerHTML = `
    <span class="node-type-badge ${typeClass}">${typeLabel}</span>
    <div class="node-name">${label || id}</div>
    ${smiles ? `<div class="smiles-box" title="Click to select all">${smiles}</div>` : ''}
    ${molblock ? `<div class="mol-3d-viewer" id="${viewerId}"></div>` : ''}
  `;

  if (molblock) {
    const viewerDiv = document.getElementById(viewerId);
    if (viewerDiv) {
      create3DViewer(viewerDiv, molblock, null);
    }
  } else if (smiles) {
    // Fallback: 2D SVG from RDKit
    const svg = await smilesTo2DSVG(smiles, 220, 150);
    if (svg) {
      const preview = document.createElement('div');
      preview.className = 'node-mol-preview';
      preview.innerHTML = svg;
      const nameEl = sidebar.querySelector('.node-name');
      if (nameEl) nameEl.insertAdjacentElement('afterend', preview);
    }
  }
}

/**
 * Render an overlay sidebar showing multiple selected nodes in one 3Dmol viewer.
 * Source = blue, target = green, intermediate = amber.
 */
function renderOverlaySidebar(sidebarId, nodesData) {
  const sidebar = document.getElementById(sidebarId);
  if (!sidebar) return;

  const hasMolblocks = nodesData.some(d => d.molblock);
  const viewerId = `sidebar-overlay-${Date.now()}`;

  const namesList = nodesData.map(d => {
    const cls = d.is_source ? 'source' : d.is_target ? 'target' : 'intermediate';
    return `<span class="node-type-badge ${cls}" style="margin:2px">${d.label || d.id}</span>`;
  }).join('');

  sidebar.innerHTML = `
    <div class="node-type-badge intermediate">${nodesData.length} nodes selected</div>
    <div class="overlay-names">${namesList}</div>
    ${hasMolblocks ? `<div class="mol-3d-viewer" id="${viewerId}"></div>` : ''}
  `;

  if (hasMolblocks && typeof $3Dmol !== 'undefined') {
    const viewerDiv = document.getElementById(viewerId);
    if (viewerDiv) {
      try {
        viewerDiv.style.width = '100%';
        viewerDiv.style.height = '220px';
        viewerDiv.style.position = 'relative';
        const viewer = $3Dmol.createViewer(viewerDiv, { backgroundColor: 'white' });
        nodesData.forEach((d, i) => {
          if (!d.molblock) return;
          viewer.addModel(d.molblock, 'sdf');
          const color = d.is_source ? '#0969da' : d.is_target ? '#1a7f37' : '#9a6700';
          viewer.setStyle({ model: i }, { stick: { color }, sphere: { scale: 0.2, color } });
        });
        viewer.zoomTo();
        viewer.render();
      } catch (e) {
        console.warn('3Dmol overlay error:', e);
      }
    }
  }
}

/* ── About modal ───────────────────────────────────────────────────────────── */
document.getElementById('about-btn').onclick = () =>
  document.getElementById('about-modal').classList.add('visible');

document.getElementById('about-close').onclick = () =>
  document.getElementById('about-modal').classList.remove('visible');

document.getElementById('about-modal').addEventListener('click', (e) => {
  if (e.target === document.getElementById('about-modal'))
    document.getElementById('about-modal').classList.remove('visible');
});

document.querySelectorAll('.modal-tab').forEach(tab => {
  tab.addEventListener('click', () => {
    const key = tab.dataset.about;
    document.querySelectorAll('.modal-tab').forEach(t => t.classList.remove('active'));
    document.querySelectorAll('.modal-section').forEach(s => s.classList.remove('active'));
    tab.classList.add('active');
    document.getElementById(`about-${key}`).classList.add('active');
  });
});

document.getElementById('copy-bibtex-btn').addEventListener('click', () => {
  const text = document.getElementById('bibtex-text').textContent;
  const btn = document.getElementById('copy-bibtex-btn');
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(() => {
      btn.textContent = 'Copied!';
      setTimeout(() => { btn.textContent = 'Copy'; }, 2000);
    }).catch(() => {
      fallbackCopy(text, btn);
    });
  } else {
    fallbackCopy(text, btn);
  }
});

/* ── Session ID + fetch wrapper ────────────────────────────────────────────── */
function getSessionId() {
  const MAX_AGE_MS = 90 * 24 * 60 * 60 * 1000; // 90 days
  let sid = localStorage.getItem('pairmap_session_id');
  const ts  = Number(localStorage.getItem('pairmap_session_ts')) || 0;
  if (!sid || (Date.now() - ts > MAX_AGE_MS)) {
    sid = crypto.randomUUID();
    localStorage.setItem('pairmap_session_id', sid);
    localStorage.setItem('pairmap_session_ts', String(Date.now()));
  }
  return sid;
}

function apiFetch(url, options = {}) {
  const headers = new Headers(options.headers || {});
  headers.set('X-Session-Id', getSessionId());
  return fetch(url, { ...options, headers });
}

function fallbackCopy(text, btn) {
  const ta = document.createElement('textarea');
  ta.value = text;
  ta.style.position = 'fixed';
  ta.style.opacity = '0';
  document.body.appendChild(ta);
  ta.select();
  try {
    document.execCommand('copy');
    btn.textContent = 'Copied!';
    setTimeout(() => { btn.textContent = 'Copy'; }, 2000);
  } catch { /* silent */ }
  document.body.removeChild(ta);
}

/**
 * Render edge MCS info in the sidebar.
 * Fetches /api/pair/{sessionId}/mcs/{nodeA}/{nodeB} and displays SVGs + legend.
 */
async function renderEdgeSidebar(sidebarId, nodeA, nodeB, similarity, sessionId) {
  const sidebar = document.getElementById(sidebarId);
  if (!sidebar || !sessionId) return;

  sidebar.innerHTML = `
    <div class="sidebar-placeholder">
      <div class="spinner" style="margin:auto"></div>
      <span style="font-size:12px;color:var(--text-muted);margin-top:0.5rem;">Computing MCS…</span>
    </div>`;

  try {
    const res = await fetch(API_BASE + `/api/pair/${sessionId}/mcs/${nodeA}/${nodeB}`);
    if (!res.ok) throw new Error(await res.text());
    const data = await res.json();

    const labelA = data.label_a || `Node ${nodeA}`;
    const labelB = data.label_b || `Node ${nodeB}`;
    const simText = similarity != null ? `Similarity: ${Number(similarity).toFixed(3)}` : '';
    sidebar.innerHTML = `
      <div class="edge-mcs-info">
        <div class="node-name" style="font-size:13px;">Edge: ${labelA} → ${labelB}</div>
        ${simText ? `<div style="font-size:12px;color:var(--text-muted);">${simText}</div>` : ''}
        <div class="mcs-legend">
          <span class="mcs-legend-item mcs-common">■ Common (${data.n_common})</span>
          <span class="mcs-legend-item mcs-deleted">■ Deleted (${data.n_deleted})</span>
          <span class="mcs-legend-item mcs-inserted">■ Inserted (${data.n_inserted})</span>
        </div>
        <div class="mcs-mol-pair">
          <div class="mcs-mol-item">
            <div class="mcs-mol-label">${labelA}</div>
            <div class="mcs-mol-svg">${data.svg_a}</div>
          </div>
          <div class="mcs-mol-item">
            <div class="mcs-mol-label">${labelB}</div>
            <div class="mcs-mol-svg">${data.svg_b}</div>
          </div>
        </div>
      </div>`;
  } catch (e) {
    sidebar.innerHTML = `
      <div class="sidebar-placeholder" style="color:var(--red);">
        MCS computation failed:<br><small>${e.message}</small>
      </div>`;
  }
}

