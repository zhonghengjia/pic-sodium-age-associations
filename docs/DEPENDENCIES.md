# Dependencies

Recorded analysis environment: R 4.6.1; mice 3.19.0, mitml 0.4-5, digest 0.6.39, jsonlite 2.0.0. Final R figure rendering recorded ggplot2 4.0.3, patchwork 1.3.2 and ragg 1.5.2. Other scripts use survival, svglite and base/grid packages. Arial was used for figures; Cairo PDF and ragg PNG devices are referenced in source.

Python code requires Python 3.11+ for hashlib.file_digest, numpy and scipy for the linear-programming separation screen, and PyMuPDF/pymupdf, Pillow and numpy for the bounded vector-text cleanup. Dependency versions are reported as recorded evidence, not installed or upgraded by this preparation.

Original scripts invoke PowerShell (`pwsh`) and local helper entry points for subprocesses. Set `PIC_PYTHON_HELPER` and `PIC_R_HELPER` to suitable local helpers when using the original guarded workflow. Optional `PIC_R_LIBRARY` selects an existing local R library. State/calibration layout QA also references external helpers supplied via `PIC_PANEL_ALIGNMENT_R` and `PIC_PANEL_ALIGNMENT_AUDIT`; those third-party/local helper implementations are not bundled. Their absence prevents the original layout-QA step from completing and must not be called a pass.

No dependency installation, lockfile generation, scientific script execution or environment migration was performed for this archive.
