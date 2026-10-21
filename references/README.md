# Saved reference inputs

`reference_snapshot.tar.gz` contains 1,904 saved API, GMT and eQTL input files. `snapshot_manifest.json` records the original paths, file sizes and SHA-256 values. `env/reference_snapshot.sha256` pins the archive itself. Open Targets responses use release 26.09.

The complete submission ZIP includes this archive. In a Git checkout, `scripts/restore_references.py` downloads the identical archive from the v1.4-submission release. The download script runs it automatically. Restoration checks the archive and every file, rejects path traversal and will not overwrite input files that differ. Network access is needed only if the archive is missing.

Primary raw datasets stay as public downloads and are not included in this snapshot. See the source/accession list and the download script.
