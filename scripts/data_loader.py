"""
data_loader.py
==============
Carga y limpieza de datos morfométricos de Atelognathus reverberii.
Compartido entre scripts de análisis y módulos de la app.
"""

import pandas as pd
import numpy as np
from pathlib import Path


SEXO_MAP = {
    'macho': 'M', 'male': 'M',
    'hembra': 'F', 'female': 'F', 'hemba': 'F',
    'juvenil': 'J', 'juvuenil': 'J', 'juvenile': 'J',
    'metamorfo': 'J',
}

CEI_DATA = {
    'reference': 'Cei (1969, 1980)',
    'males': {
        'n': 5,
        'LHU_mean': 36.6,
        'LHU_min': 35.0,
        'LHU_max': 38.0,
        'note': 'SVL (snout-vent to cloaca). Holotype + 4 paratypes.'
    },
    'females': {
        'n': 2,
        'LHU_mean': 37.25,
        'LHU_min': 36.5,
        'LHU_max': 38.0,
        'note': 'SVL. Allotype + 1 paratype. Cei 1980 reports max 45 mm.'
    },
    'note_methodology': (
        'Cei used snout-vent length (SVL, measured to cloaca). '
        'Present study uses snout-urostyle length (LHU, measured to urostyle '
        'along dorsal midline). LHU is systematically shorter than SVL; '
        'direct comparisons should be interpreted with this caveat.'
    )
}


def load_morpho(filepath: str | Path) -> pd.DataFrame:
    """
    Load morphometric data from either:
      - The clean CSV shipped with this repository
        (``data/morphometrics_laguna_azul.csv``)
      - The original Excel file (``data/Planilla_CMR.xlsx``)

    Returns a standardised DataFrame with columns:
        ID, Fecha/date, Peso/mass_g, LHU/SUL_mm, AB/MW_mm,
        Sexo/sex_class, tiene_contenido_estomacal/stomach_content_present
    All morphometric columns are float; Sexo is M/F/J or NaN.
    Only rows with complete measurements (LHU, Peso, AB) are returned.
    """
    filepath = Path(filepath)
    suffix = filepath.suffix.lower()

    if suffix == '.csv':
        df = pd.read_csv(filepath)
        # Normalise column names: accept both English and Spanish headers
        col_map = {
            'sul_mm': 'LHU', 'SUL_mm': 'LHU',
            'mass_g': 'Peso', 'Mass_g': 'Peso',
            'mw_mm': 'AB', 'MW_mm': 'AB',
            'sex_class': 'Sexo', 'Sex_class': 'Sexo',
            'stomach_content_present': 'tiene_contenido_estomacal',
            'date': 'Fecha',
        }
        df = df.rename(columns=col_map)
        # sex_class CSV uses full words; remap to single-letter codes
        if 'Sexo' in df.columns:
            df['Sexo'] = df['Sexo'].map(
                lambda x: SEXO_MAP.get(str(x).lower().strip(), np.nan)
                if pd.notna(x) else np.nan
            )
        if 'tiene_contenido_estomacal' not in df.columns:
            df['tiene_contenido_estomacal'] = False

    elif suffix in ('.xlsx', '.xls'):
        df = pd.read_excel(filepath,
                           sheet_name='Primera Captura (Marcado)', header=0)
        df.columns = [
            'Fecha', 'Hora', 'ID', 'Marca', 'Peso', 'LHU', 'AB',
            'Lavado', 'Cant_Ep', 'Sexo', 'Obs',
            'MAI', 'MAD', 'MPI', 'MPD',
            'T1', 'T2', 'T3', 'Notas'
        ]
        df = df[['ID', 'Fecha', 'Peso', 'LHU', 'AB', 'Sexo', 'Cant_Ep', 'Obs']].copy()
        for col in ['Peso', 'LHU', 'AB']:
            df[col] = pd.to_numeric(
                df[col].astype(str).str.replace(',', '.', regex=False),
                errors='coerce'
            )
        raw = df['Sexo'].astype(str).str.lower().str.strip()
        df['Sexo'] = raw.map(lambda x: SEXO_MAP.get(x, np.nan))

        def _has_content(val):
            s = str(val).strip().lower()
            if s in ('nan', '', 'negativo', 'negative', '0'):
                return False
            try:
                return float(s) > 0
            except ValueError:
                return True
        df['tiene_contenido_estomacal'] = df['Cant_Ep'].apply(_has_content)
    else:
        raise ValueError(f"Unsupported file type: {suffix}. Use .csv or .xlsx")

    df = df.dropna(subset=['LHU', 'Peso', 'AB']).reset_index(drop=True)
    return df


def load_recapturas(filepath: str | Path) -> pd.DataFrame:
    """
    Carga la hoja 'Recapturas medidas' con todos los eventos de recaptura.
    """
    df = pd.read_excel(
        filepath,
        sheet_name='Recapturas medidas',
        header=0
    )
    df.columns = [
        'ID', 'Marca', 'MAI', 'MAD', 'MPI', 'MPD',
        'Fecha', 'Hora', 'Peso', 'LHU', 'AB', 'Lavado', 'Cant_Ep',
        'Sexo', 'Obs'
    ] + [f'extra_{i}' for i in range(df.shape[1] - 15)]

    for col in ['Peso', 'LHU', 'AB']:
        df[col] = pd.to_numeric(
            df[col].astype(str).str.replace(',', '.', regex=False),
            errors='coerce'
        )
    return df


def load_cmr_matrix(filepath: str | Path) -> pd.DataFrame:
    """
    Carga la hoja 'Historial de capturauras con fe' como matriz binaria
    ID × fecha para modelos de captura-recaptura.
    """
    df = pd.read_excel(
        filepath,
        sheet_name='Historial de capturauras con fe',
        header=1
    )
    df = df.rename(columns={df.columns[0]: 'ID'})
    df = df.dropna(subset=['ID'])
    df['ID'] = df['ID'].astype(int)
    return df
