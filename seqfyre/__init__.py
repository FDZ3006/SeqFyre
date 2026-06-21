"""
SeqFyre — pipeline analisis struktur data bioinformatika.

Paket inti berorientasi objek:
    SequenceRecord : model satu sekuens (frekuensi nukleotida, GC content).
    Parser         : pembaca FASTA / FASTQ / ZIP.
    Analyzer       : pengurutan GC, statistik, CSV, dan grafik EDA.
    build_result_zip : pengemas hasil (CSV + grafik) menjadi ZIP in-memory.
"""

from .models import SequenceRecord
from .parser import Parser, ParserError
from .analyzer import Analyzer
from .packager import build_result_zip

__all__ = [
    "SequenceRecord",
    "Parser",
    "ParserError",
    "Analyzer",
    "build_result_zip",
]

__version__ = "1.0.0"
