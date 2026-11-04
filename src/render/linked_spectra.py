"""The precursor (MS1) spectrum of the selected spectrum.

A fragmentation spectrum's precursor spectrum is the MS1 survey spectrum its
precursor was isolated from. The precursor spectrum cell shows that spectrum
next to the selected one and highlights the precursor mass in it, so the
precursor's charge states and isotope envelope can be checked by eye. For an
MS1 spectrum it shows the spectrum itself.

The link comes from FLASHDeconv's results: `PrecursorScan` of the 3D S/N plot
data names the precursor spectrum, else the last MS1 spectrum before the
fragmentation spectrum is taken. The precursor mass is the mass of the
precursor spectrum within 0.01 Da of the fragmentation spectrum's precursor
mass, as the 3D S/N plot pairs them. Like insight_data, the module is free of
Streamlit and openms-insight.
"""

import polars as pl

from src.render import insight_data

# A mass of the precursor spectrum within this of the precursor mass is the precursor
PRECURSOR_MASS_TOLERANCE_DA = 1e-2

LINK_SCHEMA = {
    'scan_idx': pl.Int64,
    'survey_idx': pl.Int64,
    'survey_mass_idx': pl.Int64,
}


def link_surveys(scans, precursor_scans=None, survey_masses=None, tolerance=PRECURSOR_MASS_TOLERANCE_DA):
    """One row per spectrum: its index (`scan_idx`), the index of its
    precursor spectrum (`survey_idx`, itself for an MS1 spectrum, null when
    there is none) and the index of the precursor mass among that spectrum's
    masses (`survey_mass_idx`, null when FLASHDeconv did not find it there).

    scans: `index`, `Scan`, `MSLevel` and `PrecursorMass` per spectrum.
    precursor_scans: `index` and the scan number of its precursor spectrum
        (`PrecursorScan`, 0 or less when unknown); optional.
    survey_masses: `scan_idx`, `mass_idx` and `MonoMass` per deconvolved mass;
        optional.
    """
    if scans.is_empty():
        return pl.DataFrame(schema=LINK_SCHEMA)
    scans = (
        scans.select(
            pl.col('index').cast(pl.Int64).alias('scan_idx'),
            pl.col('Scan').cast(pl.Int64),
            pl.col('MSLevel').cast(pl.Int64),
            pl.col('PrecursorMass').cast(pl.Float64),
        )
        .sort('scan_idx')
    )
    # The last MS1 spectrum at or before each spectrum
    scans = scans.with_columns(
        pl.when(pl.col('MSLevel') == 1).then(pl.col('scan_idx')).forward_fill().alias('previous_ms1')
    )
    if precursor_scans is not None and not precursor_scans.is_empty():
        named = (
            precursor_scans.select(
                pl.col('index').cast(pl.Int64).alias('scan_idx'),
                pl.col('PrecursorScan').cast(pl.Int64),
            )
            .filter(pl.col('PrecursorScan') > 0)
            .unique('scan_idx', keep='first')
            .join(
                scans.select(pl.col('Scan').alias('PrecursorScan'), pl.col('scan_idx').alias('named_idx'))
                .unique('PrecursorScan', keep='first'),
                on='PrecursorScan',
            )
            .select('scan_idx', 'named_idx')
        )
        scans = scans.join(named, on='scan_idx', how='left')
    else:
        scans = scans.with_columns(pl.lit(None, dtype=pl.Int64).alias('named_idx'))
    links = scans.select(
        'scan_idx',
        pl.when(pl.col('MSLevel') == 1).then(pl.col('scan_idx'))
        .otherwise(pl.coalesce('named_idx', 'previous_ms1')).alias('survey_idx'),
        pl.when(pl.col('MSLevel') > 1).then(pl.col('PrecursorMass')).alias('PrecursorMass'),
    )
    if survey_masses is None or survey_masses.is_empty():
        return links.select(
            'scan_idx', 'survey_idx', pl.lit(None, dtype=pl.Int64).alias('survey_mass_idx')
        )
    # The precursor mass among the masses of the precursor spectrum: the closest
    # within the tolerance
    found = (
        links.filter(pl.col('PrecursorMass') > 0)
        .join(
            survey_masses.select(
                pl.col('scan_idx').cast(pl.Int64).alias('survey_idx'),
                pl.col('mass_idx').cast(pl.Int64).alias('survey_mass_idx'),
                pl.col('MonoMass').cast(pl.Float64),
            ),
            on='survey_idx',
        )
        .with_columns((pl.col('MonoMass') - pl.col('PrecursorMass')).abs().alias('error'))
        .filter(pl.col('error') < tolerance)
        .sort(['scan_idx', 'error', 'survey_mass_idx'])
        .unique('scan_idx', keep='first')
        .select('scan_idx', 'survey_mass_idx')
    )
    return links.join(found, on='scan_idx', how='left').select(list(LINK_SCHEMA))


def survey_links(file_manager, dataset_id):
    """link_surveys() of a dataset."""
    scans = insight_data.scan_table(file_manager, dataset_id).collect()
    precursor_scans = None
    if insight_data.result_path(file_manager, dataset_id, 'threedim_SN_plot') is not None:
        threedim = insight_data.scan(file_manager, dataset_id, 'threedim_SN_plot')
        if 'PrecursorScan' in threedim.collect_schema().names():
            precursor_scans = threedim.select('index', 'PrecursorScan').collect()
    survey_masses = None
    if insight_data.result_path(file_manager, dataset_id, 'mass_table') is not None:
        survey_masses = (
            insight_data.mass_table(file_manager, dataset_id)
            .select('scan_idx', 'mass_idx', 'MonoMass')
            .collect()
        )
    return link_surveys(scans, precursor_scans, survey_masses)


def survey_of(links, scan_idx):
    """(survey_idx, survey_mass_idx) of a spectrum, or (None, None)."""
    if scan_idx is None:
        return None, None
    row = links.filter(pl.col('scan_idx') == int(scan_idx))
    if row.is_empty():
        return None, None
    return row['survey_idx'][0], row['survey_mass_idx'][0]
