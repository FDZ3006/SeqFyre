"""
parser.py
---------
Kelas Parser: membaca berkas FASTA, FASTQ, dan ZIP
lalu mengubahnya menjadi List berisi objek SequenceRecord.

Strategi parsing berlapis:
    1. Jika BioPython terpasang -> pakai Bio.SeqIO
    2. Jika tidak               -> pakai parser manual pure-Python
"""

import io
import zipfile

from .models import SequenceRecord

# Cek ketersediaan BioPython saat modul dimuat
try:
    from Bio import SeqIO
    _HAS_BIOPYTHON = True
except ImportError:
    _HAS_BIOPYTHON = False


class ParserError(Exception):
    # Dilempar ketika berkas tidak bisa dibaca atau formatnya tidak dikenal
    pass


class Parser:
    # Ekstensi berkas yang didukung
    FASTA_EXT = (".fasta", ".fa", ".fna", ".ffn", ".frn", ".txt")
    FASTQ_EXT = (".fastq", ".fq")

    def __init__(self, use_biopython=True):
        # Bisa dipaksa nonaktif untuk pengujian fallback
        self.use_biopython = use_biopython and _HAS_BIOPYTHON

    def get_engine(self):
        # Nama mesin parsing yang sedang dipakai
        if self.use_biopython:
            return "BioPython"
        return "Fallback (pure-Python)"

    # Tambahkan property engine agar kompatibel dengan kode lama
    engine = property(get_engine)

    def parse_path(self, path):
        # Parse berkas dari path di disk, mendukung .fasta/.fastq/.zip
        lower = path.lower()
        if lower.endswith(".zip"):
            with open(path, "rb") as fh:
                return self.parse_zip(fh.read())
        fmt = self._format_from_name(lower)
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return self._parse_text(fh.read(), fmt)

    def parse_upload(self, filename, data):
        # Parse berkas hasil upload Flask berupa bytes dan nama berkas
        lower = filename.lower()
        if lower.endswith(".zip"):
            return self.parse_zip(data)
        fmt = self._format_from_name(lower)
        text = data.decode("utf-8", errors="replace")
        return self._parse_text(text, fmt)

    def parse_zip(self, data):
        # Ekstrak semua berkas FASTA/FASTQ di dalam ZIP lalu gabungkan
        records = []
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            members = []
            for m in zf.namelist():
                # Lewati direktori dan file sampah macOS
                if m.endswith("/"):
                    continue
                if "__MACOSX" in m:
                    continue
                if m.lower().endswith(self.FASTA_EXT + self.FASTQ_EXT):
                    members.append(m)
            if not members:
                raise ParserError("ZIP tidak memuat berkas FASTA/FASTQ yang dikenali.")
            for member in members:
                fmt = self._format_from_name(member.lower())
                content = zf.read(member).decode("utf-8", errors="replace")
                records.extend(self._parse_text(content, fmt))
        return records

    def _format_from_name(self, lower_name):
        # Tentukan format dari ekstensi berkas
        if lower_name.endswith(self.FASTQ_EXT):
            return "fastq"
        if lower_name.endswith(self.FASTA_EXT):
            return "fasta"
        # Default: anggap FASTA
        return "fasta"

    def _parse_text(self, text, fmt):
        # Pilih mesin parsing lalu kembalikan daftar SequenceRecord
        text = text.strip()
        if not text:
            raise ParserError("Berkas kosong - tidak ada sekuens untuk dibaca.")

        # Auto-koreksi format jika tidak sesuai ekstensi
        if text.startswith("@") and fmt != "fastq":
            fmt = "fastq"
        elif text.startswith(">") and fmt != "fasta":
            fmt = "fasta"

        if self.use_biopython:
            records = self._parse_biopython(text, fmt)
        else:
            records = self._parse_fallback(text, fmt)

        if not records:
            raise ParserError("Tidak ada sekuens valid yang berhasil dibaca dari berkas.")
        return records

    def _parse_biopython(self, text, fmt):
        # Parsing menggunakan Bio.SeqIO
        handle = io.StringIO(text)
        out = []
        for rec in SeqIO.parse(handle, fmt):
            description = rec.description or rec.id
            out.append(SequenceRecord(rec.id, description, str(rec.seq)))
        return out

    def _parse_fallback(self, text, fmt):
        # Parser manual tanpa dependensi eksternal
        if fmt == "fastq":
            return self._parse_fastq_manual(text)
        return self._parse_fasta_manual(text)

    def _parse_fasta_manual(self, text):
        # Baca FASTA baris per baris, gabungkan sekuens multi-baris
        records = []
        header = None
        seq_parts = []

        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Simpan rekaman sebelumnya
                if header is not None:
                    records.append(self._build(header, seq_parts))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)

        # Simpan rekaman terakhir
        if header is not None:
            records.append(self._build(header, seq_parts))
        return records

    def _parse_fastq_manual(self, text):
        # Baca FASTQ: 4 baris per rekaman (@header, seq, +, quality)
        records = []
        lines = [ln for ln in text.splitlines() if ln.strip() != ""]
        i = 0
        while i < len(lines):
            if not lines[i].startswith("@"):
                i += 1
                continue
            if i + 1 >= len(lines):
                break
            header = lines[i][1:].strip()
            seq = lines[i + 1].strip()
            records.append(self._build(header, [seq]))
            i += 4  # lewati baris + dan quality
        return records

    def _build(self, header, seq_parts):
        # Bentuk SequenceRecord dari header dan bagian-bagian sekuens
        seq = "".join(seq_parts)
        seq_id = header.split()[0] if header else "unknown"
        return SequenceRecord(seq_id, header, seq)
