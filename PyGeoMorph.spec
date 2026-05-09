# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller specification for PyGeoMorph.
Builds a single-file executable that bundles the Flask app, the HTML
template, and the scientific Python stack.

NOTE on macOS: we deliberately build for the current native architecture
only (no universal2). When run on a macos-13 GitHub runner this produces
an x86_64 binary; on macos-14+ runners (Apple Silicon) it produces an
arm64 binary. Building universal2 fails because numpy/scipy/opencv wheels
on PyPI are normally single-arch.
"""

import sys
from pathlib import Path
from PyInstaller.utils.hooks import collect_submodules, collect_data_files

block_cipher = None
HERE = Path(SPECPATH).resolve() if 'SPECPATH' in globals() else Path('.').resolve()

# ---------------------------------------------------------------------
# Hidden imports — modules that PyInstaller's static analysis misses.
# ---------------------------------------------------------------------
hiddenimports = []
# numpy: collect tutti i submoduli per evitare i problemi di hook
hiddenimports += collect_submodules('numpy')
# scipy: i submoduli più usati
hiddenimports += collect_submodules('scipy')
hiddenimports += collect_submodules('scipy.special')
hiddenimports += collect_submodules('scipy.linalg')
hiddenimports += collect_submodules('scipy.spatial')
hiddenimports += collect_submodules('scipy.stats')
# Flask + estensioni
hiddenimports += collect_submodules('flask')
hiddenimports += [
    'cv2',
    'PIL.Image',
    'PIL.ImageFilter',
    'pdf2image',
    'jinja2.ext',
    'werkzeug.middleware.proxy_fix',
    # numpy internals che alcuni hook mancano
    'numpy.core._methods',
    'numpy.core._dtype_ctypes',
    'numpy.lib.format',
]

# ---------------------------------------------------------------------
# Bundled data files: the Flask templates folder.
# ---------------------------------------------------------------------
datas = [
    ('templates', 'templates'),
]
try:
    datas += collect_data_files('numpy')
except Exception:
    pass
try:
    datas += collect_data_files('scipy')
except Exception:
    pass

# ---------------------------------------------------------------------
a = Analysis(
    ['app.py'],
    pathex=[str(HERE)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'tkinter', 'matplotlib', 'IPython', 'jupyter',
        'pytest', 'sphinx', 'pyqt5', 'pyqt6', 'PySide2', 'PySide6',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# Determine the icon (optional)
icon_path = None
for candidate in ('icon.ico', 'icon.icns', 'icon.png'):
    if (HERE / candidate).exists():
        icon_path = str(HERE / candidate)
        break

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipped_data,
    a.datas,
    [],
    name='PyGeoMorph',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=icon_path,
)
