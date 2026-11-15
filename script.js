"use strict";

const CONTRIBUTIONS = {
  "-CH3": 23.58, "-CH2-": 22.88, ">CH-": 21.74,
  "-OH": 92.88, "-O-": 22.42, "-COOH": 169.06,
  "CH aromatique": 26.73, "-Cl": 38.13
};

const BASE_MOLECULES = {
  "CCO":      { nom: "éthanol",        groupes: { "-CH3": 1, "-CH2-": 1, "-OH": 1 } },
  "CCCO":     { nom: "propanol",       groupes: { "-CH3": 1, "-CH2-": 2, "-OH": 1 } },
  "c1ccccc1": { nom: "benzène",        groupes: { "CH aromatique": 6 } },
  "CC(=O)O":  { nom: "acide acétique", groupes: { "-CH3": 1, "-COOH": 1 } }
};

const $ = id => document.getElementById(id);
const count = (s, ch) => s.split(ch).length - 1;

/* ---------- Conversion en °C ---------- */
function toCelsius(str) {
  if (!str) return [null, null];
  const s = String(str);
  const m = s.match(/(-?\d+\.?\d*)/);
  if (!m) return [null, s];
  const val = parseFloat(m[1]);
  if (isNaN(val)) return [null, s];
  if (s.includes("°F") || s.includes(" F at") || (s.includes("F ") && !s.includes("°C"))) {
    const c = (val - 32) * 5 / 9;
    return [c, c.toFixed(1) + " °C"];
  }
  if (s.includes(" K") || s.includes("°K") || (val > 200 && !s.includes("°C") && !s.includes("°F"))) {
    if (val > 150) { const c = val - 273.15; return [c, c.toFixed(1) + " °C"]; }
    return [val, s];
  }
  return [val, val.toFixed(1) + " °C"];
}

/* ---------- Estimation de secours (Clapeyron / Trouton) ---------- */
function estimationSecours(smiles) {
  const s = smiles.trim();
  if (["[Fe]", "Fe", "[fe]"].includes(s)) return [3134.0, "Clapeyron-Métal (Fe)"];
  if (["[He]", "He"].includes(s)) return [4.2, "Clapeyron-Gaz rare (He)"];
  if (["N#N", "N2", "[N]#N", "N=N"].includes(s)) return [77.3, "Clapeyron-Diatomique (N2)"];
  if (["O", "[O]", "O=O", "O2"].includes(s)) return [373.15, "Clapeyron-Diatomique (O2/H2O)"];
  if (["[Ar]", "Ar"].includes(s)) return [87.3, "Clapeyron-Gaz rare (Ar)"];

  const nbC = count(s, "C") + count(s, "c");
  const nbO = count(s, "O");
  const nbN = count(s, "N");
  const nbCl = count(s, "Cl");
  if (nbC + nbO + nbN + nbCl > 0) {
    return [198.2 + nbC * 21.5 + nbO * 45 + nbN * 18 + nbCl * 30, "Corrélation étendue (Clapeyron)"];
  }
  return [198.2 + s.length * 20, "Estimation brute"];
}

/* ---------- PubChem ---------- */
const PUG = "https://pubchem.ncbi.nlm.nih.gov/rest";

async function fetchJSON(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error("HTTP " + r.status);
  return r.json();
}
async function fetchText(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error("HTTP " + r.status);
  return r.text();
}

function scanBoiling(obj, depth = 0) {
  if (depth > 25 || obj === null) return null;
  if (Array.isArray(obj)) {
    for (const it of obj) { const r = scanBoiling(it, depth + 1); if (r) return r; }
    return null;
  }
  if (typeof obj === "object") {
    if (String(obj.TOCHeading || "").includes("Boiling")) {
      for (const info of obj.Information || []) {
        const v = info.Value || {};
        for (const key of ["StringWithMarkup", "String", "Number"]) {
          if (key in v) {
            const val = v[key];
            if (Array.isArray(val) && val.length) {
              return typeof val[0] === "object" ? (val[0].String || "") : String(val[0]);
            }
            return String(val);
          }
        }
      }
    }
    for (const k of Object.keys(obj)) { const r = scanBoiling(obj[k], depth + 1); if (r) return r; }
  }
  return null;
}

async function obtenirInfosPubChem(smiles) {
  const corr = { Fe: "[Fe]", He: "[He]", N2: "N#N", O2: "O=O" };
  const sm = corr[smiles] || smiles;
  try {
    const cidData = await fetchJSON(`${PUG}/pug/compound/smiles/${encodeURIComponent(sm)}/cids/JSON`);
    const cid = cidData.IdentifierList.CID[0];

    let tbTxt = null, tbNum = null;
    try {
      const data = await fetchJSON(`${PUG}/pug_view/data/compound/${cid}/JSON`);
      let raw = scanBoiling(data);
      if (!raw) {
        const m = JSON.stringify(data).match(/(-?\d+\.?\d*\s*°[CF]\s*at\s*760\s*mmHg)/);
        if (m) raw = m[1];
      }
      if (raw) [tbNum, tbTxt] = toCelsius(raw);
    } catch (e) { console.warn("Point d'ébullition:", e); }
    if (!tbTxt) tbTxt = `PAS DE DONNÉE pour CID ${cid}`;

    let sdf = null;
    try { sdf = await fetchText(`${PUG}/pug/compound/cid/${cid}/SDF?record_type=3d`); }
    catch { try { sdf = await fetchText(`${PUG}/pug/compound/cid/${cid}/SDF`); } catch { sdf = null; } }

    return { cid, tbTxt, tbNum, sdf, png: `${PUG}/pug/compound/cid/${cid}/PNG?image_size=600x600` };
  } catch (e) {
    console.warn("PubChem:", e);
    return null;
  }
}

/* ---------- Lecture du SDF ---------- */
function parserSDF(text) {
  if (!text) return { atomes: [], liaisons: [] };
  const L = text.split(/\r?\n/);
  let idx = L.findIndex(l => l.includes("V2000"));
  if (idx < 0) idx = 3;
  const nA = parseInt(L[idx].substr(0, 3)), nL = parseInt(L[idx].substr(3, 3));
  if (isNaN(nA) || isNaN(nL)) return { atomes: [], liaisons: [] };
  const atomes = [], liaisons = [];
  for (let i = idx + 1; i < idx + 1 + nA; i++) {
    const x = parseFloat(L[i].substr(0, 10)), y = parseFloat(L[i].substr(10, 10)), z = parseFloat(L[i].substr(20, 10));
    const sym = L[i].substr(31, 3).trim();
    if (![x, y, z].some(isNaN)) atomes.push({ x, y, z, sym });
  }
  for (let i = idx + 1 + nA; i < idx + 1 + nA + nL; i++) {
    const a = parseInt(L[i].substr(0, 3)) - 1, b = parseInt(L[i].substr(3, 3)) - 1, o = parseInt(L[i].substr(6, 3));
    if (!isNaN(a) && !isNaN(b)) liaisons.push({ a, b, o });
  }
  return { atomes, liaisons };
}

/* ---------- Visualiseur 3D (canvas, sans bibliothèque) ---------- */
const COLORS = { C: "#2b2b2b", O: "#d62828", N: "#2f5bea", H: "#cfd6d8", Cl: "#2e9e4f", S: "#e0c200", Fe: "#e07b00", He: "#19b5c9" };
const RADII = { H: 6, He: 9 };
const viewer = { atomes: [], liaisons: [], yaw: 0.6, pitch: 0.4, zoom: 1, scale: 1, raf: 0 };

function setupViewer(atomes, liaisons) {
  const cx = atomes.reduce((s, a) => s + a.x, 0) / atomes.length;
  const cy = atomes.reduce((s, a) => s + a.y, 0) / atomes.length;
  const cz = atomes.reduce((s, a) => s + a.z, 0) / atomes.length;
  viewer.atomes = atomes.map(a => ({ ...a, x: a.x - cx, y: a.y - cy, z: a.z - cz }));
  viewer.liaisons = liaisons;
  const maxR = Math.max(1, ...viewer.atomes.map(a => Math.hypot(a.x, a.y, a.z)));
  viewer.scale = 190 / maxR;
  viewer.yaw = 0.6; viewer.pitch = 0.4; viewer.zoom = 1;
  draw();
}

function draw() {
  const cv = $("viewer"), ctx = cv.getContext("2d");
  const W = cv.width, H = cv.height;
  ctx.clearRect(0, 0, W, H);
  const cyw = Math.cos(viewer.yaw), syw = Math.sin(viewer.yaw);
  const cp = Math.cos(viewer.pitch), sp = Math.sin(viewer.pitch);
  const k = viewer.scale * viewer.zoom;
  const P = viewer.atomes.map(a => {
    const x1 = a.x * cyw + a.z * syw, z1 = -a.x * syw + a.z * cyw;
    const y2 = a.y * cp - z1 * sp, z2 = a.y * sp + z1 * cp;
    return { X: W / 2 + x1 * k, Y: H / 2 - y2 * k, Z: z2, sym: a.sym };
  });
  const items = [];
  viewer.liaisons.forEach(l => {
    const A = P[l.a], B = P[l.b];
    if (A && B) items.push({ t: "b", z: (A.Z + B.Z) / 2, A, B, o: l.o });
  });
  P.forEach(p => items.push({ t: "a", z: p.Z, p }));
  items.sort((u, v) => u.z - v.z);
  ctx.lineCap = "round";
  for (const it of items) {
    if (it.t === "b") {
      ctx.strokeStyle = "#8a9a9e";
      ctx.lineWidth = it.o > 1 ? 7 : 4;
      ctx.beginPath(); ctx.moveTo(it.A.X, it.A.Y); ctx.lineTo(it.B.X, it.B.Y); ctx.stroke();
      if (it.o > 1) {
        ctx.strokeStyle = "#e6ecec"; ctx.lineWidth = 2.5; ctx.stroke();
      }
    } else {
      const s = it.p.sym, r = (RADII[s] || 13) * Math.min(1.6, viewer.zoom);
      ctx.fillStyle = COLORS[s] || "#e07b00";
      ctx.beginPath(); ctx.arc(it.p.X, it.p.Y, r, 0, Math.PI * 2); ctx.fill();
      ctx.strokeStyle = "#111"; ctx.lineWidth = 1; ctx.stroke();
      if (s !== "H") {
        ctx.fillStyle = getComputedStyle(document.body).color;
        ctx.font = "bold 13px sans-serif";
        ctx.fillText(s, it.p.X + r + 2, it.p.Y - r);
      }
    }
  }
}

(function bindViewer() {
  const cv = $("viewer");
  let drag = null;
  cv.addEventListener("pointerdown", e => { drag = { x: e.clientX, y: e.clientY }; cv.setPointerCapture(e.pointerId); });
  cv.addEventListener("pointerup", () => { drag = null; });
  cv.addEventListener("pointermove", e => {
    if (!drag) return;
    viewer.yaw += (e.clientX - drag.x) * 0.01;
    viewer.pitch += (e.clientY - drag.y) * 0.01;
    drag = { x: e.clientX, y: e.clientY };
    draw();
  });
  cv.addEventListener("wheel", e => {
    e.preventDefault();
    viewer.zoom = Math.min(4, Math.max(0.3, viewer.zoom * (e.deltaY < 0 ? 1.1 : 0.9)));
    draw();
  }, { passive: false });
})();

/* ---------- Analyse ---------- */
function setStatus(msg, err = false) {
  const s = $("status");
  s.textContent = msg;
  s.className = err ? "err" : "";
}

async function analyser(smiles) {
  $("go").disabled = true;
  $("result").hidden = true;
  setStatus("Recherche en cours sur PubChem…");

  let nom = smiles, methode = "Joback", tbK, groupes = null, total = 0;
  if (BASE_MOLECULES[smiles]) {
    const info = BASE_MOLECULES[smiles];
    nom = info.nom; groupes = info.groupes;
    for (const [g, n] of Object.entries(groupes)) total += (CONTRIBUTIONS[g] || 0) * n;
    tbK = 198.2 + total;
  } else {
    [tbK, methode] = estimationSecours(smiles);
  }
  const tbC = tbK - 273.15;

  const info = await obtenirInfosPubChem(smiles);

  $("nom").textContent = nom;
  $("cid").textContent = info ? `${smiles} · CID ${info.cid}` : `${smiles} · introuvable sur PubChem`;
  $("methode").textContent = methode;
  $("tbEst").textContent = `${tbC.toFixed(1)} °C (${tbK.toFixed(1)} K)`;
  $("tbReel").textContent = info ? info.tbTxt : "Non disponible";
  let err = "—";
  if (info && info.tbNum !== null && info.tbNum !== 0) {
    err = (Math.abs(tbC - info.tbNum) / Math.abs(info.tbNum) * 100).toFixed(1) + " %";
  }
  $("erreur").textContent = err;

  const tab = $("tabJoback");
  if (groupes) {
    tab.innerHTML = "<tr><th>Groupe</th><th>Calcul</th><th>Contribution (K)</th></tr>" +
      Object.entries(groupes).map(([g, n]) =>
        `<tr><td>${g}</td><td>${n} × ${CONTRIBUTIONS[g]}</td><td>${(CONTRIBUTIONS[g] * n).toFixed(2)}</td></tr>`).join("") +
      `<tr><td colspan="2">Tb = 198,2 + somme</td><td>${tbK.toFixed(2)}</td></tr>`;
    $("detailJoback").hidden = false;
  } else {
    $("detailJoback").hidden = true;
  }

  const img = $("img2d");
  if (info) { img.src = info.png; img.hidden = false; } else { img.hidden = true; }

  const { atomes, liaisons } = parserSDF(info && info.sdf);
  const cv = $("viewer");
  if (atomes.length) {
    cv.parentElement.hidden = false;
    setupViewer(atomes, liaisons);
  } else {
    cv.parentElement.hidden = true;
  }

  $("result").hidden = false;
  setStatus(info ? "" : "Molécule non trouvée sur PubChem : seule l'estimation est affichée.", !info);
  $("go").disabled = false;
}

$("form").addEventListener("submit", e => {
  e.preventDefault();
  const s = $("smiles").value.trim();
  if (s) analyser(s);
});
$("chips").addEventListener("click", e => {
  const s = e.target.dataset && e.target.dataset.s;
  if (s) { $("smiles").value = s; analyser(s); }
});
