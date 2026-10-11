"""
Returned-soil reference compositions for the Apollo/Luna landing sites
(average soils; Lunar Sourcebook Table 7.15; Basu & Riegsecker 1998),
read from data/reference/validation.xlsx.
"""
import openpyxl

# Landing-site coordinates (deg, lon in [-180, 180])
SITE_COORDS = {
    "11": (0.674, 23.473),
    "12": (-3.012, -23.42),
    "14": (-3.645, -17.471),
    "15": (26.132, 3.633),
    "16": (8.973, 15.499),
    "17": (20.19, 30.77),
    "L16": (-0.68, 56.3),
    "L20": (3.53, 56.55),
    "L24": (12.714, 62.2),
}

ELEMENT_ROWS = {"Ti": 33, "Al": 34, "Fe": 35, "Mg": 36, "Ca": 37, "Si": 38}  # openpyxl 1-indexed row
SITE_COLS = {  # 0-indexed column in the sheet
    "11": 1, "12": 2, "14": 3, "15a": 4, "15b": 5, "15c": 6, "15": 7,
    "16a": 8, "16b": 9, "16c": 10, "16": 11, "17a": 12, "17b": 13, "17c": 14,
    "17d": 15, "17": 16, "L16": 17, "L20": 18, "L24": 19,
}


def load_reference(xlsx_path):
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb["Sheet1"]
    ref = {}
    for site, coords in SITE_COORDS.items():
        col = SITE_COLS[site]
        wt = {el: ws.cell(row=r, column=col + 1).value for el, r in ELEMENT_ROWS.items()}
        if wt.get("Si") in (None, 0):
            continue
        ratios = {el: wt[el] / wt["Si"] for el in wt if el != "Si"}
        ref[site] = {"lat": coords[0], "lon": coords[1], "wt_pct": wt, "ratio_to_Si": ratios}
    return ref
