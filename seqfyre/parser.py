"""
parser.py
=========
Kelas Parser: membaca berkas FASTA / FASTQ (dan arsip ZIP berisi keduanya)
lalu mengubahnya menjadi List[SequenceRecord].

Strategi parsing berlapis:
    1. Jika BioPython terpasang  -> pakai Bio.SeqIO (sesuai ekspektasi bonus).
    2. Jika tidak                -> pakai parser manual pure-Python (fallback).

Dengan begitu aplikasi tetap berjalan di lingkungan minim dependensi,
namun tetap memanfaatkan BioPython saat tersedia (mis. di Render).
"""

from __future__ import annotations

import io
import zipfile
from typing import Iterable

from .models import SequenceRecord

# Deteksi ketersediaan BioPython sekali saat modul dimuat.
try:
    from Bio import SeqIO  # type: ignore

    _HAS_BIOPYTHON = True
except ImportError:  # pragma: no cover - tergantung lingkungan
    _HAS_BIOPYTHON = False


class ParserError(Exception):
    """Dilempar ketika berkas tidak dapat di-parse atau formatnya tak dikenal."""


class Parser:
    """Mengubah teks/berkas FASTA-FASTQ-ZIP menjadi List[SequenceRecord]."""

    FASTA_EXT = (".fasta", ".fa", ".fna", ".ffn", ".frn", ".txt")
    FASTQ_EXT = (".fastq", ".fq")

    def __init__(self, use_biopython: bool = True) -> None:
        # Boleh dipaksa nonaktif (mis. untuk pengujian fallback).
        self.use_biopython = use_biopython and _HAS_BIOPYTHON

    @property
    def engine(self) -> str:
        """Nama mesin parsing yang sedang dipakai (untuk pelaporan)."""
        return "BioPython" if self.use_biopython else "Fallback (pure-Python)"

    # ------------------------------------------------------------------ #
    # API publik
    # ------------------------------------------------------------------ #
    def parse_path(self, path: str) -> list[SequenceRecord]:
        """Parse berkas dari path pada disk. Mendukung .fasta/.fastq/.zip."""
        lower = path.lower()
        if lower.endswith(".zip"):
            with open(path, "rb") as fh:
                return self.parse_zip(fh.read())
        fmt = self._format_from_name(lower)
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return self._parse_text(fh.read(), fmt)

    def parse_upload(self, filename: str, data: bytes) -> list[SequenceRecord]:
        """Parse berkas hasil upload (Flask) berupa bytes + nama berkas asli."""
        lower = filename.lower()
        if lower.endswith(".zip"):
            return self.parse_zip(data)
        fmt = self._format_from_name(lower)
        text = data.decode("utf-8", errors="replace")
        return self._parse_text(text, fmt)

    def parse_zip(self, data: bytes) -> list[SequenceRecord]:
        """Ekstrak semua berkas FASTA/FASTQ di dalam ZIP lalu gabungkan."""
        records: list[SequenceRecord] = []
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            members = [
                m for m in zf.namelist()
                if not m.endswith("/")                 # lewati direktori
                and "__MACOSX" not in m                 # lewati sampah macOS
                and m.lower().endswith(self.FASTA_EXT + self.FASTQ_EXT)
            ]
            if not members:
                raise ParserError(
                    "ZIP tidak memuat berkas FASTA/FASTQ yang dikenali."
                )
            for member in members:
                fmt = self._format_from_name(member.lower())
                content = zf.read(member).decode("utf-8", errors="replace")
                records.extend(self._parse_text(content, fmt))
        return records

    # ------------------------------------------------------------------ #
    # Helper internal
    # ------------------------------------------------------------------ #
    def _format_from_name(self, lower_name: str) -> str:
        """Tentukan 'fasta' atau 'fastq' dari ekstensi berkas."""
        if lower_name.endswith(self.FASTQ_EXT):
            return "fastq"
        if lower_name.endswith(self.FASTA_EXT):
            return "fasta"
        # Default aman: anggap FASTA (format paling umum & toleran).
        return "fasta"

    def _parse_text(self, text: str, fmt: str) -> list[SequenceRecord]:
        """Pilih mesin parsing lalu kembalikan daftar SequenceRecord."""
        text = text.strip()
        if not text:
            raise ParserError("Berkas kosong — tidak ada sekuens untuk dibaca.")

        # Auto-koreksi: kadang berkas .fastq sebetulnya FASTA, atau sebaliknya.
        if text.startswith("@") and fmt != "fastq":
            fmt = "fastq"
        elif text.startswith(">") and fmt != "fasta":
            fmt = "fasta"

        if self.use_biopython:
            records = self._parse_biopython(text, fmt)
        else:
            records = self._parse_fallback(text, fmt)

        if not records:
            raise ParserError(
                "Tidak ada sekuens valid yang berhasil dibaca dari berkas."
            )
        return records

    def _parse_biopython(self, text: str, fmt: str) -> list[SequenceRecord]:
        """Parsing memakai Bio.SeqIO."""
        handle = io.StringIO(text)
        out: list[SequenceRecord] = []
        for rec in SeqIO.parse(handle, fmt):
            description = rec.description or rec.id
            out.append(SequenceRecord(rec.id, description, str(rec.seq)))
        return out

    def _parse_fallback(self, text: str, fmt: str) -> list[SequenceRecord]:
        """Parser manual tanpa dependensi eksternal."""
        if fmt == "fastq":
            return self._parse_fastq_manual(text)
        return self._parse_fasta_manual(text)

    @staticmethod
    def _parse_fasta_manual(text: str) -> list[SequenceRecord]:
        """Baca FASTA baris-demi-baris; gabungkan sekuens multi-baris."""
        records: list[SequenceRecord] = []
        header: str | None = None
        seq_parts: list[str] = []

        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Simpan rekaman sebelumnya sebelum mulai yang baru.
                if header is not None:
                    records.append(Parser._build(header, seq_parts))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
        # Jangan lupa rekaman terakhir.
        if header is not None:
            records.append(Parser._build(header, seq_parts))
        return records

    @staticmethod
    def _parse_fastq_manual(text: str) -> list[SequenceRecord]:
        """Baca FASTQ (4 baris per rekaman: @header, seq, +, quality)."""
        records: list[SequenceRecord] = []
        lines = [ln for ln in text.splitlines() if ln.strip() != ""]
        i = 0
        while i + 3 < len(lines) or (i + 1 < len(lines) and lines[i].startswith("@")):
            if not lines[i].startswith("@"):
                i += 1
                continue
            header = lines[i][1:].strip()
            seq = lines[i + 1].strip()
            records.append(Parser._build(header, [seq]))
            i += 4  # lompat ke rekaman berikutnya (lewati '+' dan quality)
        return records

    @staticmethod
    def _build(header: str, seq_parts: Iterable[str]) -> SequenceRecord:
        """Bentuk SequenceRecord; kata pertama header = accession/SeqID."""
        seq = "".join(seq_parts)
        seq_id = header.split()[0] if header else "unknown"
        return SequenceRecord(seq_id, header, seq)
