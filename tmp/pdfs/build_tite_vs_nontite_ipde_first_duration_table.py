from pathlib import Path

import build_tite_vs_nontite_duration_table as base


ROOT = Path(__file__).resolve().parents[2]
PRESENTATION = ROOT / "Presentation 8-17-2026"
RAW_DATA = PRESENTATION / "Raw Data"

base.OUTPUT_DIR = PRESENTATION / "Table and Plots" / "TITE"
base.OUTPUT_PDF = base.OUTPUT_DIR / (
    "phase12_tite_vs_nontite_3designs_N30_fut0p85_"
    "scenarios_1_16_20_24_27_38_tables_ipde_first.pdf"
)
base.TITE_FILE = RAW_DATA / "TITE_AIDE_phase_I_II_IDX_1001_to_2000_dose_summary.csv"
base.TITE_ONE_STAGE_REPLACEMENT_FILE = RAW_DATA / (
    "TITE_AIDE_phase_I_II_IDX_1001_to_2000_newdesign_dose_summary.csv"
)
base.NON_TITE_FILE = RAW_DATA / (
    "AIDE_phase_I_II_modelsadditive_shared_N30_ncycle2_"
    "rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_"
    "ap0p15x0p85_IDX_0001_to_1000_dose_summary_ipdefirst.csv"
)
base.NON_TITE_ONE_STAGE_REPLACEMENT_FILE = RAW_DATA / (
    "AIDE_phase_I_II_modelsadditive_shared_N30_ncycle2_"
    "rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_"
    "ap0p15x0p85_IDX_0001_to_1000_newdesign_ipdefirst_dose_summary.csv"
)
base.NON_TITE_LABEL = "IPDE-first non-TITE"
base.SUBTITLE = (
    "TITE versus IPDE-first non-TITE: one-stage, two-stage highest utility, "
    "and two-stage top-2 randomized"
)
base.SOURCE_NOTE = (
    "Sources: updated TITE and IPDE-first non-TITE one-stage newdesign summaries; "
    "existing two-stage summaries."
)
base.SIMULATION_NOTE = (
    "Both sources use 1,000 simulations per populated row; "
    "TITE No OBD % = 100 - sum of dose OBD selections."
)
base.DURATION_NOTE = "Duration is mean calendar days to complete the trial."


if __name__ == "__main__":
    base.main()
