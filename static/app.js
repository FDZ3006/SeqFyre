/* ============================================================
   SeqFyre — frontend logic
   Menangani: unggah berkas, pilihan Top-N, panggilan ke backend,
   render hasil (statistik, kolom termal, tabel, grafik), unduhan,
   dan pengalih bahasa ID/EN.
   ============================================================ */

// -------------------- Kamus terjemahan (i18n) -------------------- //
const I18N = {
  id: {
    tagline: "16S rRNA · analisis termal",
    eyebrow: "Pipeline struktur data bioinformatika",
    title: 'Urutkan sekuens.<br /><span class="accent">Temukan panasnya.</span>',
    subtitle:
      "Unggah berkas FASTA/FASTQ gen 16S rRNA. SeqFyre menghitung GC content, " +
      "mengurutkannya, dan menyingkap gradien termal — bakteri termofil naik ke " +
      "puncak, mesofil mengendap di dasar.",
    dz_main: "Letakkan berkas di sini atau klik untuk memilih",
    dz_sub: "FASTA · FASTQ · ZIP — maks 16 MB",
    topn_label: "Tampilkan teratas:",
    demo_btn: "Coba dataset demo",
    analyze_btn: "Analisis",
    download_btn: "Unduh ZIP (CSV + grafik)",
    sec_summary: "Ringkasan",
    sec_thermal: "Kolom termal — sekuens terbaik",
    sec_table: "Tabel lengkap (terurut GC)",
    sec_eda: "Eksplorasi data (EDA)",
    foot: "SeqFyre · Mini Project Struktur Data Bioinformatika (BIF1223) · IPB University",
    // status & dinamis
    busy: "Menganalisis…",
    done: "Selesai — mesin: ",
    no_file: "Pilih berkas terlebih dahulu.",
    file_ready: "Berkas siap: ",
    err: "Galat: ",
    // label statistik
    st_total: "Total sekuens",
    st_meangc: "Rata-rata GC",
    st_thermo: "Termofil",
    st_meso: "Mesofil",
    // header tabel
    th_rank: "Rank", th_acc: "Aksesi", th_org: "Organisme", th_len: "Panjang",
    th_gc: "GC (%)", th_a: "A", th_t: "T", th_g: "G", th_c: "C", th_n: "N", th_cls: "Kelas",
    cls_thermo: "Termofil", cls_meso: "Mesofil",
    // judul grafik
    p_hist: "Histogram distribusi GC",
    p_bar: "GC per sekuens",
    p_comp: "Komposisi nukleotida",
    p_scatter: "Panjang vs GC",
  },
  en: {
    tagline: "16S rRNA · thermal analysis",
    eyebrow: "Bioinformatics data-structure pipeline",
    title: 'Sort the sequences.<br /><span class="accent">Find the heat.</span>',
    subtitle:
      "Upload a 16S rRNA FASTA/FASTQ file. SeqFyre computes GC content, sorts it, " +
      "and reveals the thermal gradient — thermophiles rise to the top, mesophiles " +
      "settle at the bottom.",
    dz_main: "Drop a file here or click to choose",
    dz_sub: "FASTA · FASTQ · ZIP — max 16 MB",
    topn_label: "Show top:",
    demo_btn: "Try demo dataset",
    analyze_btn: "Analyze",
    download_btn: "Download ZIP (CSV + charts)",
    sec_summary: "Summary",
    sec_thermal: "Thermal column — best sequences",
    sec_table: "Full table (sorted by GC)",
    sec_eda: "Exploratory data analysis (EDA)",
    foot: "SeqFyre · Bioinformatics Data Structures Mini Project (BIF1223) · IPB University",
    busy: "Analyzing…",
    done: "Done — engine: ",
    no_file: "Choose a file first.",
    file_ready: "File ready: ",
    err: "Error: ",
    st_total: "Total sequences",
    st_meangc: "Mean GC",
    st_thermo: "Thermophiles",
    st_meso: "Mesophiles",
    th_rank: "Rank", th_acc: "Accession", th_org: "Organism", th_len: "Length",
    th_gc: "GC (%)", th_a: "A", th_t: "T", th_g: "G", th_c: "C", th_n: "N", th_cls: "Class",
    cls_thermo: "Thermophile", cls_meso: "Mesophile",
    p_hist: "GC distribution histogram",
    p_bar: "GC per sequence",
    p_comp: "Nucleotide composition",
    p_scatter: "Length vs GC",
  },
};

let lang = "id";
let lastResult = null;     // payload terakhir (untuk render ulang saat ganti bahasa)
let selectedFile = null;
let topN = 3;

const $ = (sel) => document.querySelector(sel);
const t = (key) => I18N[lang][key] ?? key;

// -------------------- i18n: terapkan ke DOM -------------------- //
function applyI18n() {
  document.documentElement.lang = lang;
  document.querySelectorAll("[data-i18n]").forEach((el) => {
    el.innerHTML = t(el.dataset.i18n);
  });
  // Render ulang hasil agar label tabel/grafik ikut berganti bahasa.
  if (lastResult) renderResults(lastResult);
}

document.querySelectorAll(".lang button").forEach((btn) => {
  btn.addEventListener("click", () => {
    lang = btn.dataset.lang;
    document.querySelectorAll(".lang button").forEach((b) => b.classList.remove("active"));
    btn.classList.add("active");
    applyI18n();
  });
});

// -------------------- Pemilihan berkas -------------------- //
const dropzone = $("#dropzone");
const fileInput = $("#fileInput");

dropzone.addEventListener("click", () => fileInput.click());
dropzone.addEventListener("keydown", (e) => {
  if (e.key === "Enter" || e.key === " ") { e.preventDefault(); fileInput.click(); }
});
fileInput.addEventListener("change", () => setFile(fileInput.files[0]));

["dragover", "dragenter"].forEach((ev) =>
  dropzone.addEventListener(ev, (e) => { e.preventDefault(); dropzone.classList.add("drag"); })
);
["dragleave", "drop"].forEach((ev) =>
  dropzone.addEventListener(ev, (e) => { e.preventDefault(); dropzone.classList.remove("drag"); })
);
dropzone.addEventListener("drop", (e) => {
  if (e.dataTransfer.files.length) setFile(e.dataTransfer.files[0]);
});

function setFile(file) {
  if (!file) return;
  selectedFile = file;
  dropzone.classList.add("has-file");
  dropzone.querySelector(".dz-main").textContent = file.name;
  $("#analyzeBtn").disabled = false;
  setStatus(t("file_ready") + file.name, "ok");
}

// -------------------- Top-N -------------------- //
$("#topnSeg").addEventListener("click", (e) => {
  const btn = e.target.closest("button");
  if (!btn) return;
  topN = parseInt(btn.dataset.n, 10);
  $("#topnSeg").querySelectorAll("button").forEach((b) => b.classList.remove("active"));
  btn.classList.add("active");
  // Jika sudah ada hasil, perbarui kolom termal mengikuti Top-N baru.
  if (lastResult && selectedFile) analyze();
  else if (lastResult && lastResult._demo) runDemo();
});

// -------------------- Status helper -------------------- //
function setStatus(msg, cls = "") {
  const el = $("#status");
  el.textContent = msg;
  el.className = "status " + cls;
}

// -------------------- Panggilan backend -------------------- //
$("#analyzeBtn").addEventListener("click", analyze);
$("#demoBtn").addEventListener("click", runDemo);

async function analyze() {
  if (!selectedFile) { setStatus(t("no_file"), "err"); return; }
  const fd = new FormData();
  fd.append("file", selectedFile);
  fd.append("top_n", topN);
  await postAnalyze("/analyze", fd, false);
}

async function runDemo() {
  const fd = new FormData();
  fd.append("top_n", topN);
  await postAnalyze("/demo", fd, true);
}

async function postAnalyze(url, formData, isDemo) {
  setStatus(t("busy"), "busy");
  $("#analyzeBtn").disabled = true;
  $("#demoBtn").disabled = true;
  try {
    const res = await fetch(url, { method: "POST", body: formData });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || res.statusText);
    data._demo = isDemo;
    lastResult = data;
    renderResults(data);
    setStatus(t("done") + data.engine, "ok");
    $("#results").classList.add("show");
  } catch (err) {
    setStatus(t("err") + err.message, "err");
  } finally {
    $("#analyzeBtn").disabled = !selectedFile;
    $("#demoBtn").disabled = false;
  }
}

// -------------------- Render hasil -------------------- //
function renderResults(data) {
  renderStats(data.summary);
  renderThermal(data.top);
  renderTable(data.table);
  renderPlots(data.plots);
  $("#downloadBtn").onclick = () => {
    window.location.href = "/download/" + data.token;
  };
}

function renderStats(s) {
  $("#stats").innerHTML = `
    <div class="stat"><div class="num">${s.total_sequences}</div><div class="cap">${t("st_total")}</div></div>
    <div class="stat"><div class="num">${s.mean_gc}%</div><div class="cap">${t("st_meangc")}</div></div>
    <div class="stat"><div class="num hot">${s.thermophile_count}</div><div class="cap">${t("st_thermo")}</div></div>
    <div class="stat"><div class="num cool">${s.mesophile_count}</div><div class="cap">${t("st_meso")}</div></div>`;
}

function renderThermal(top) {
  const max = Math.max(...top.map((r) => r["GC_Content(%)"]), 1);
  $("#thermal").innerHTML = top.map((r, i) => {
    const isThermo = r.Classification === "Termofil";
    const cls = isThermo ? "thermo" : "meso";
    const clsName = isThermo ? t("cls_thermo") : t("cls_meso");
    const width = (r["GC_Content(%)"] / max) * 100;
    return `
      <div class="band ${cls} ${i === 0 ? "top" : ""}">
        <div class="fill" style="width:${width}%"></div>
        <span class="rank">${String(i + 1).padStart(2, "0")}</span>
        <span class="acc">${r.Accession_ID}</span>
        <span class="org">${r.Organism}</span>
        <span class="pill">${clsName}</span>
        <span class="gc">${r["GC_Content(%)"]}%</span>
      </div>`;
  }).join("");
}

function renderTable(rows) {
  const heads = ["th_rank","th_acc","th_org","th_len","th_gc","th_a","th_t","th_g","th_c","th_n","th_cls"];
  const thead = `<thead><tr>${heads.map((h) => `<th>${t(h)}</th>`).join("")}</tr></thead>`;
  const tbody = "<tbody>" + rows.map((r) => {
    const isThermo = r.Classification === "Termofil";
    const clsName = isThermo ? t("cls_thermo") : t("cls_meso");
    const clsCss = isThermo ? "cls-thermo" : "cls-meso";
    return `<tr>
      <td class="mono">${r.Rank}</td>
      <td class="mono">${r.Accession_ID}</td>
      <td><i>${r.Organism}</i></td>
      <td class="mono">${r.Length}</td>
      <td class="mono">${r["GC_Content(%)"]}</td>
      <td class="mono">${r.A_Count}</td>
      <td class="mono">${r.T_Count}</td>
      <td class="mono">${r.G_Count}</td>
      <td class="mono">${r.C_Count}</td>
      <td class="mono">${r.N_Count}</td>
      <td class="mono ${clsCss}">${clsName}</td>
    </tr>`;
  }).join("") + "</tbody>";
  $("#resultTable").innerHTML = thead + tbody;
}

function renderPlots(plots) {
  const order = [
    ["histogram_gc", "p_hist"],
    ["barchart_gc", "p_bar"],
    ["composition", "p_comp"],
    ["scatter_len_gc", "p_scatter"],
  ];
  $("#plots").innerHTML = order.map(([key, label]) =>
    plots[key]
      ? `<div class="plot"><h4>${t(label)}</h4><img src="${plots[key]}" alt="${t(label)}" /></div>`
      : ""
  ).join("");
}

// init
applyI18n();
