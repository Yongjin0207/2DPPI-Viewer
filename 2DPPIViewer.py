import argparse
import json
import re
from collections import defaultdict
from html import escape


# -----------------------------
# Residue classes
# -----------------------------
POSITIVE = {"ARG", "LYS", "HIS"}
NEGATIVE = {"ASP", "GLU"}
POLAR_UNCHARGED = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
SPECIAL = {"GLY", "PRO"}
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP"}

# Only truly charged atom names
POS_CHARGE_ATOMNAMES = {
    "LYS": {"NZ"},
    "ARG": {"NE", "NH1", "NH2"},
    # If you do NOT want His highlighted, remove this entry:
    "HIS": {"ND1", "NE2"},
}
NEG_CHARGE_ATOMNAMES = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}

# (Changed) colors for +charge / -charge atoms
POS_CHARGE_COLOR = "#FF6D00"  # orange
NEG_CHARGE_COLOR = "#6A1B9A"  # purple
DISULPHIDE_COLOR = "#FFD54F"  # yellow


# -----------------------------
# UI strings
# -----------------------------
UI = {
    "title": "2DPPI-Viewer",
    "btn_drag_on": "Drag nodes: ON",
    "btn_drag_off": "Drag nodes: OFF",
    "btn_lock_on": "Chain-side lock: ON",
    "btn_lock_off": "Chain-side lock: OFF",
    "btn_pan_on": "Pan & Zoom: ON",
    "btn_pan_off": "Pan & Zoom: OFF",
    "btn_contacts_on": "Contacts: ON",
    "btn_contacts_off": "Contacts: OFF",
    "btn_hbonly_on": "H-bond residues only: ON",
    "btn_hbonly_off": "H-bond residues only: OFF",
    "btn_fit": "Fit to view",
    "btn_reset": "Reset",
    "btn_export": "Export SVG",
    "hint": (
        "Drag residues • Wheel zoom • Pan: Shift-drag or middle mouse on background • "
        "Alt + drag (focus residue) rotates its atom diagram • "
        "Green dashed = H-bond (distance) • Yellow dashed = disulphide • Gray = contact"
    ),
}


# -----------------------------
# Small helpers
# -----------------------------
def clean_res_label(line: str) -> str:
    m = re.search(r"([A-Za-z]{3}\d+\([^)]+\))\s*$", line.strip())
    if m:
        return m.group(1)
    toks = line.split()
    return toks[-1] if toks else line.strip()


def is_res_label(s: str) -> bool:
    s = s.strip()
    return bool(re.search(r"[A-Za-z]", s) and re.search(r"\d", s) and not re.fullmatch(r"\d+(\.\d+)?", s))


def slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", s)


def display_res_name(res_id: str) -> str:
    return re.sub(r"\([^)]+\)\s*$", "", res_id).strip()


def extract_chain(res_id: str) -> str:
    m = re.search(r"\(([^)]+)\)\s*$", res_id)
    if not m:
        return "?"
    ch = m.group(1).strip()
    return ch[:1] if ch else "?"


def res3_from_label(res_id: str) -> str:
    name = display_res_name(res_id)
    m = re.match(r"([A-Za-z]{3})", name)
    return (m.group(1).upper() if m else "UNK")


def atom_color(el: str) -> str:
    el = (el or "").upper()
    if el == "C":
        return "#222"
    if el == "O":
        return "#D32F2F"
    if el == "N":
        return "#1976D2"
    if el == "S":
        return "#F9A825"
    if el == "P":
        return "#8E24AA"
    return "#555"


def residue_class(res3: str) -> str:
    r = res3.upper()
    if r in POSITIVE:
        return "positive"
    if r in NEGATIVE:
        return "negative"
    if r in POLAR_UNCHARGED:
        return "polar"
    if r in SPECIAL:
        return "special"
    if r in HYDROPHOBIC:
        return "hydrophobic"
    return "other"



AA1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}

def label1_from_label(label3: str) -> str:
    """Convert 'LYS123' -> 'K123' (best-effort)."""
    m = re.match(r"^([A-Za-z]{3})(.*)$", (label3 or "").strip())
    if not m:
        return label3
    res3 = m.group(1).upper()
    rest = m.group(2)
    return AA1.get(res3, "X") + rest

def class_colors(cls: str):
    if cls == "positive":
        return ("#ff6d60", "#b23b32")
    if cls == "negative":
        return ("#6aa9ff", "#2e5aa8")
    if cls == "polar":
        return ("#ffd54f", "#c49000")
    if cls == "special":
        return ("#99F6E4", "#0F766E")  # mint fill, teal stroke
    if cls == "hydrophobic":
        return ("#d0d0d0", "#7a7a7a")
    return ("#95A5A6", "#666666")


def override_atom_color(res3: str, atom_name: str, base_color: str) -> str:
    res3 = (res3 or "").upper()
    atom_name = (atom_name or "").upper()
    if res3 in POS_CHARGE_ATOMNAMES and atom_name in POS_CHARGE_ATOMNAMES[res3]:
        return POS_CHARGE_COLOR
    if res3 in NEG_CHARGE_ATOMNAMES and atom_name in NEG_CHARGE_ATOMNAMES[res3]:
        return NEG_CHARGE_COLOR
    if res3 == 'CYS' and atom_name == 'SG':
        return '#FFD54F'
    return base_color


# -----------------------------
# DRW parser
# -----------------------------
def parse_drw(path: str):
    lines = open(path, "r", encoding="utf-8", errors="ignore").read().splitlines()

    residues = []
    atom_by_idx = {}

    i = 0
    while i < len(lines):
        if lines[i].strip().startswith("#R"):
            i += 1
            if i >= len(lines):
                break

            res_label = clean_res_label(lines[i].strip())
            if not is_res_label(res_label):
                i += 1
                continue

            chain = extract_chain(res_label)

            while i < len(lines) and lines[i].strip() != "#A":
                i += 1
            if i >= len(lines):
                break
            i += 1  # past "#A"

            atoms = []
            while i < len(lines) and not lines[i].strip().startswith("#"):
                l = lines[i]
                mid = re.search(r"\[(\d+)\]\s*$", l)
                if mid:
                    idx = int(mid.group(1))
                    parts = l.split()
                    if len(parts) >= 6:
                        try:
                            x = float(parts[1])
                            y = float(parts[2])
                        except Exception:
                            i += 1
                            continue

                        name = parts[0]     # atom name (NH1, OE1, CA, ...)
                        element = parts[5]  

                        atom_by_idx[idx] = {
                            "idx": idx,
                            "name": name,
                            "x": x,
                            "y": y,
                            "element": element,
                            "res": res_label,
                            "chain": chain,
                        }
                        atoms.append(idx)
                i += 1

            residues.append({"id": res_label, "chain": chain, "atoms": atoms})
            continue

        i += 1

    bonds, hbonds, contacts, disulphides = [], [], [], []
    for l in lines:
        m = re.match(r"^\s*([0125])\s+(\d+)\s+(\d+)(.*)$", l)
        if not m:
            continue
        t = int(m.group(1))
        a = int(m.group(2))
        b = int(m.group(3))
        rest = m.group(4).strip()

        if a not in atom_by_idx or b not in atom_by_idx:
            continue

        if t == 0:
            order = rest.split()[0] if rest else "1"
            bonds.append((a, b, order))
        elif t == 1:
            dist = None
            if rest:
                try:
                    dist = float(rest.split()[0])
                except Exception:
                    dist = None
            hbonds.append((a, b, dist))
        elif t == 2:
            contacts.append((a, b))
        elif t == 5:
            disulphides.append((a, b))

    return residues, atom_by_idx, bonds, hbonds, contacts, disulphides


# -----------------------------
# Geometry builders
# -----------------------------
def compute_residue_centers(res_atoms, atom_by_idx):
    centers = {}
    for rid, atoms in res_atoms.items():
        xs = [atom_by_idx[a]["x"] for a in atoms if a in atom_by_idx]
        ys = [atom_by_idx[a]["y"] for a in atoms if a in atom_by_idx]
        centers[rid] = (sum(xs) / len(xs), sum(ys) / len(ys)) if xs else (0.0, 0.0)
    return centers


def compute_atom_local_offsets(res_center, atom_by_idx):
    atom_local = defaultdict(dict)
    for idx, info in atom_by_idx.items():
        rid = info["res"]
        cx, cy = res_center[rid]
        atom_local[rid][idx] = (info["x"] - cx, info["y"] - cy, info["element"], info["name"])
    return atom_local


def group_covalent_bonds_by_residue(bonds, atom_by_idx):
    bonds_by_res = defaultdict(list)
    for a, b, order in bonds:
        ra = atom_by_idx[a]["res"]
        rb = atom_by_idx[b]["res"]
        if ra == rb:
            bonds_by_res[ra].append((a, b, order))
    return bonds_by_res


def collect_focus_residues_hbonds_and_disulphides(hbonds, disulphides, atom_by_idx):
    focus_res = set()
    atom_hbonds = []
    atom_disulphides = []
    for a, b, dist in hbonds:
        ra = atom_by_idx[a]["res"]
        rb = atom_by_idx[b]["res"]
        if ra == rb:
            continue
        if dist is not None:
            focus_res.add(ra)
            focus_res.add(rb)
            atom_hbonds.append((a, b, dist))
    for a, b in disulphides:
        ra = atom_by_idx[a]["res"]
        rb = atom_by_idx[b]["res"]
        if ra == rb:
            continue
        focus_res.add(ra)
        focus_res.add(rb)
        atom_disulphides.append((a, b))
    return focus_res, atom_hbonds, atom_disulphides


def dedup_contact_pairs(contacts, atom_by_idx):
    contact_pairs = defaultdict(int)
    for a, b in contacts:
        ra = atom_by_idx[a]["res"]
        rb = atom_by_idx[b]["res"]
        if ra == rb:
            continue
        key = tuple(sorted([ra, rb]))
        contact_pairs[key] += 1
    return contact_pairs


# -----------------------------
# HTML builder
# -----------------------------
def build_html(
    drw_path: str,
    out_html: str,
    width: int,
    height: int,
    margin: int,
    atom_scale: float,
    split_extra: float,
    x_expand: float,
    y_expand: float,
    relax_iters: int,
    relax_gap_px: float,
    relax_step: float,
    simple_r: float,
    focus_r: float,
    relax_on_pointerup: bool,
    boundary_extend_px: float,
    atom_label_mode: str,
    show_atom_labels_for_all: bool,
    chain_label_offset_px: float,
):
    residues, atom_by_idx, bonds, hbonds, contacts, disulphides = parse_drw(drw_path)

    res_atoms = {r["id"]: r["atoms"] for r in residues}
    res_chain = {r["id"]: r["chain"] for r in residues}

    res_center = compute_residue_centers(res_atoms, atom_by_idx)
    atom_local = compute_atom_local_offsets(res_center, atom_by_idx)
    bonds_by_res = group_covalent_bonds_by_residue(bonds, atom_by_idx)

    focus_res, atom_hbonds, atom_disulphides = collect_focus_residues_hbonds_and_disulphides(hbonds, disulphides, atom_by_idx)
    contact_pairs = dedup_contact_pairs(contacts, atom_by_idx)

    # Split chains in Y
    y_vals = [y for _, y in res_center.values()]
    y_range = (max(y_vals) - min(y_vals)) if y_vals else 1.0
    split = (y_range / 2.0) + split_extra

    adj_center = {}
    for rid, (x, y) in res_center.items():
        ch = res_chain.get(rid, "?")
        y2 = y - split if ch == "A" else (y + split if ch == "B" else y)
        adj_center[rid] = (x, y2)

    yA = [adj_center[r][1] for r in adj_center if res_chain.get(r) == "A"]
    yB = [adj_center[r][1] for r in adj_center if res_chain.get(r) == "B"]
    boundary = (max(yA) + min(yB)) / 2.0 if yA and yB else 0.0

    # Bounds include atom diagram extents for focus residues
    max_local = 0.0
    for rid in focus_res:
        for _, (dx, dy, *_rest) in atom_local[rid].items():
            max_local = max(max_local, abs(dx), abs(dy))
    ext_drw = max_local * atom_scale

    minx = min(x - (ext_drw if rid in focus_res else 0) for rid, (x, _y) in adj_center.items())
    maxx = max(x + (ext_drw if rid in focus_res else 0) for rid, (x, _y) in adj_center.items())
    miny = min(y - (ext_drw if rid in focus_res else 0) for rid, (_x, y) in adj_center.items())
    maxy = max(y + (ext_drw if rid in focus_res else 0) for rid, (_x, y) in adj_center.items())

    x_span = (maxx - minx) if (maxx - minx) != 0 else 1.0
    y_span = (maxy - miny) if (maxy - miny) != 0 else 1.0
    scale = min((width - 2 * margin) / x_span, (height - 2 * margin) / y_span)

    def to_screen(x, y):
        return margin + (x - minx) * scale, margin + (y - miny) * scale

    boundary_y = to_screen(minx, boundary)[1]

    # Residue data (screen coords)
    res_data = {}
    for rid, (x, y) in adj_center.items():
        sx, sy = to_screen(x, y)
        ch = res_chain.get(rid, "?")
        focus = rid in focus_res

        label = display_res_name(rid)
        res3 = res3_from_label(rid)
        cls = residue_class(res3)

        entry = {
            "x": sx,
            "y": sy,
            "chain": ch,
            "focus": focus,
            "label3": label,
            "label1": label1_from_label(label),
            "label": label,  # initial (3-letter); JS can toggle
            "res3": res3,
            "cls": cls,
        }

        if focus:
            offs = {}
            for idx, (dx, dy, el, name) in atom_local[rid].items():
                offs[str(idx)] = {
                    "dx": dx * atom_scale * scale,
                    "dy": dy * atom_scale * scale,
                    "el": el,
                    "name": name,
                }
            entry["atoms"] = offs
            entry["bonds"] = [[str(a), str(b), order] for a, b, order in bonds_by_res[rid]]

        res_data[rid] = entry

    # Expand layout
    screen_xs = [r["x"] for r in res_data.values()]
    screen_ys = [r["y"] for r in res_data.values()]
    cx = sum(screen_xs) / len(screen_xs) if screen_xs else width / 2
    cy = sum(screen_ys) / len(screen_ys) if screen_ys else height / 2

    for r in res_data.values():
        r["x"] = cx + (r["x"] - cx) * x_expand
        r["y"] = cy + (r["y"] - cy) * y_expand

    # ---------------------------------
    # X-axis zoning: H-bond (focus) residues closer to center-x,
    # other residues farther from center-x (y is preserved).
    # ---------------------------------
    HB_INNER_FACTOR = 0.7
    OTHER_OUTER_FACTOR = 1.08
    SHELL_GAP_PX = 25.0

    def apply_layer_layout_yonly(res_data, boundary_y):
        """
        Two-layer separation around the horizontal boundary (x-axis concept):
        - focus (atom-rendered / H-bond network) residues are pulled closer to boundary_y
        - other residues are pushed farther away from boundary_y
        Only Y is modified; X is preserved.
        """
        for _rid, r in res_data.items():
            y = r["y"]
            # Determine which side of the boundary the residue is on
            side = -1.0 if y < boundary_y else 1.0  # -1: upper, +1: lower
            dist = abs(y - boundary_y)

            if r.get("focus"):
                new_dist = max(8.0, dist * HB_INNER_FACTOR)  # keep very close but not zero
            else:
                new_dist = (dist * OTHER_OUTER_FACTOR) + SHELL_GAP_PX

            r["y"] = boundary_y + side * new_dist

    apply_layer_layout_yonly(res_data, boundary_y)



    # Edge list
    edges = []
    for (u, v), cnt in contact_pairs.items():
        edges.append({"type": "contact", "u": u, "v": v, "count": cnt})
    for a, b, dist in atom_hbonds:
        ra = atom_by_idx[a]["res"]
        rb = atom_by_idx[b]["res"]
        edges.append({"type": "hbond", "ra": ra, "rb": rb, "a": str(a), "b": str(b), "dist": dist})
    for a, b in atom_disulphides:
        ra = atom_by_idx[a]["res"]
        rb = atom_by_idx[b]["res"]
        edges.append({"type": "disulphide", "ra": ra, "rb": rb, "a": str(a), "b": str(b)})

    # -----------------------------
    # SVG
    # -----------------------------
    svg_parts = []

    bx1 = -boundary_extend_px
    bx2 = width + boundary_extend_px
    svg_parts.append(
        f'<line id="boundary" x1="{bx1:.2f}" y1="{boundary_y:.2f}" x2="{bx2:.2f}" y2="{boundary_y:.2f}" '
        f'stroke="#222" stroke-dasharray="7,6" stroke-width="1.5" opacity="0.45"/>'
    )

    svg_parts.append(
        f'<text id="chainA_label" x="-200" y="{boundary_y - chain_label_offset_px:.2f}" '
        f'font-size="18" font-weight="900" fill="#333" text-anchor="middle">Chain A</text>'
    )
    svg_parts.append(
        f'<text id="chainB_label" x="-200" y="{boundary_y + chain_label_offset_px + 18:.2f}" '
        f'font-size="18" font-weight="900" fill="#333" text-anchor="middle">Chain B</text>'
    )

    # Contacts
    for e in edges:
        if e["type"] == "contact":
            eid = f'c_{slug(e["u"])}__{slug(e["v"])}'
            svg_parts.append(
                f'<line class="edge contact" id="{eid}" x1="0" y1="0" x2="0" y2="0" '
                f'stroke="#9aa0a6" stroke-width="1.2" opacity="0.55"/>'
            )

    # Hbonds / disulphides
    for e in edges:
        if e["type"] == "hbond":
            eid = f'h_{slug(e["ra"])}_{e["a"]}__{slug(e["rb"])}_{e["b"]}'
            svg_parts.append(
                f'<line class="edge hbond" id="{eid}" x1="0" y1="0" x2="0" y2="0" '
                f'stroke="#1b8e5a" stroke-width="2.4" stroke-dasharray="6,4" opacity="0.95"/>'
            )
            svg_parts.append(
                f'<text class="edgeLabel" id="{eid}_t" x="0" y="0" font-size="14" fill="#1b8e5a" '
                f'text-anchor="middle" dominant-baseline="central"></text>'
            )
        elif e["type"] == "disulphide":
            eid = f'd_{slug(e["ra"])}_{e["a"]}__{slug(e["rb"])}_{e["b"]}'
            svg_parts.append(
                f'<line class="edge disulphide" id="{eid}" x1="0" y1="0" x2="0" y2="0" '
                f'stroke="{DISULPHIDE_COLOR}" stroke-width="2.6" stroke-dasharray="6,4" opacity="0.98"/>'
            )

    # Nodes
    for rid, r in res_data.items():
        gid = f'res_{slug(rid)}'
        focus = bool(r["focus"])
        label = r["label"]
        cls = r["cls"]
        res3 = r["res3"]
        fill, stroke = class_colors(cls)
        x, y = r["x"], r["y"]
        chain = r["chain"]

        svg_parts.append(
            f'<g class="res {"focus" if focus else "simple"}" id="{gid}" '
            f'data-res="{escape(rid)}" data-chain="{escape(chain)}" data-cls="{escape(cls)}" data-res3="{escape(res3)}" '
            f'transform="translate({x:.2f},{y:.2f})">'
        )

        if focus:
            svg_parts.append(f'<g class="atomBlock" id="atomblock_{slug(rid)}">')

            # Bonds
            for a, b, order in r.get("bonds", []):
                da = r["atoms"][a]
                db = r["atoms"][b]
                x1_, y1_, x2_, y2_ = da["dx"], da["dy"], db["dx"], db["dy"]
                width2 = 2.6 if order != "1" else 2.0
                svg_parts.append(
                    f'<line class="bond" x1="{x1_:.2f}" y1="{y1_:.2f}" x2="{x2_:.2f}" y2="{y2_:.2f}" '
                    f'stroke="#333" stroke-width="{width2}" stroke-linecap="round" opacity="0.95"/>'
                )

            # Atoms (± charge override by atom name)
            for idx, ainfo in r["atoms"].items():
                ax, ay = ainfo["dx"], ainfo["dy"]
                el = (ainfo.get("el") or "").upper()
                name = (ainfo.get("name") or "").upper()

                base = atom_color(el)
                col = override_atom_color(res3, name, base)

                svg_parts.append(
                    f'<circle class="atom" data-atom="{escape(idx)}" data-el="{escape(el)}" data-name="{escape(name)}" '
                    f'cx="{ax:.2f}" cy="{ay:.2f}" r="5.3" fill="{col}" stroke="#fff" stroke-width="1.3"/>'
                )

                txt = name if (atom_label_mode == "atomname") else el
                if show_atom_labels_for_all or el in ["O", "N", "S", "P"]:
                    t2 = txt[:4] if len(txt) > 4 else txt
                    svg_parts.append(
                        f'<text class="atomLabel" x="{ax+8:.2f}" y="{ay-6:.2f}" font-size="11" fill="{col}">{escape(t2)}</text>'
                    )

            svg_parts.append(f'<circle class="center" cx="0" cy="0" r="7.0" fill="{fill}" opacity="0.22"/>')
            svg_parts.append("</g>")  # atomBlock

            svg_parts.append(
                f'<text class="resLabel" x="0" y="-40" font-size="18" font-weight="900" fill="{stroke}" '
                f'text-anchor="middle">{escape(label)}</text>'
            )
        else:
            svg_parts.append(
                f'<circle class="node" cx="0" cy="0" r="{simple_r:.1f}" fill="{fill}" fill-opacity="0.22" stroke="{stroke}" stroke-width="2"/>'
            )
            svg_parts.append(
                f'<text class="nodeLabel" x="0" y="0" font-size="12" fill="{stroke}" '
                f'text-anchor="middle" dominant-baseline="central">{escape(label)}</text>'
            )

        svg_parts.append("</g>")

    svg_str = "\n".join(svg_parts)

    # Legend (residues on left, atom-level on right)
    # Colors: use class fill with alpha for backgrounds; use class stroke for text to match node labels
    def rgba_from_hex(hex_color: str, alpha: float):
        h = hex_color.lstrip("#")
        r = int(h[0:2], 16); g = int(h[2:4], 16); b = int(h[4:6], 16)
        return f"rgba({r},{g},{b},{alpha})"

    # residue legend pills (label, fill_hex, stroke_hex)
    residue_items = [
        ("positive", "Positive", *class_colors("positive")),
        ("negative", "Negative", *class_colors("negative")),
        ("polar", "Polar", *class_colors("polar")),
        ("special", "Special", *class_colors("special")),
        ("hydrophobic", "Hydrophobic", *class_colors("hydrophobic")),
    ]
    # atom legend dots (label, color_hex)
    atom_items = [
        ("posAtom", "+charge atoms", POS_CHARGE_COLOR),
        ("negAtom", "-charge atoms", NEG_CHARGE_COLOR),
        ("disAtom", "disulphide", DISULPHIDE_COLOR),
    ]

    residue_html = "".join(
        f'<span class="lgPill" data-cls="{cls}" style="background:{rgba_from_hex(fill, 0.22)}; border-color:{stroke}; color:{stroke};">{escape(title)}</span>'
        for (cls, title, fill, stroke) in residue_items
    )
    atom_html = "".join(
        f'<span class="lgAtomItem"><span class="lgDot" data-atomkey="{k}" style="background:{c};"></span>'
        f'<span class="lgAtomText">{escape(t)}</span></span>'
        for (k, t, c) in atom_items
    )

    legend_html = (
        f'<div class="legendGroup legendResidues">{residue_html}</div>'
        f'<div class="legendSpacer"></div>'
        f'<div class="legendGroup legendAtoms">{atom_html}</div>'
    )


    # -----------------------------
    # HTML
    # -----------------------------
    js_extra = """
// -----------------------------
// Theme + AA label toggle (added)
// -----------------------------
const btnAA = document.getElementById("btnAA");
const btnTheme = document.getElementById("btnTheme");
const themePanel = document.getElementById("themePanel");
const btnThemeReset = document.getElementById("btnThemeReset");
const btnThemeClose = document.getElementById("btnThemeClose");

function $(id) { return document.getElementById(id); }

const THEME_DEFAULT = {
  positive: {fill: "#ff6d60", stroke: "#b23b32"},
  negative: {fill: "#6aa9ff", stroke: "#2e5aa8"},
  polar: {fill: "#ffd54f", stroke: "#c49000"},
  special: {fill: "#99F6E4", stroke: "#0F766E"},
  hydrophobic: {fill: "#d0d0d0", stroke: "#7a7a7a"},
  other: {fill: "#95A5A6", stroke: "#666666"},
  posAtom: "#FF6D00",
  negAtom: "#6A1B9A",
  disAtom: "#FFD54F",
  alpha: 0.22
};

function loadTheme() {
  try {
    const s = localStorage.getItem("dimplot_theme_v1");
    if (!s) return JSON.parse(JSON.stringify(THEME_DEFAULT));
    const t = JSON.parse(s);
    const out = JSON.parse(JSON.stringify(THEME_DEFAULT));
    for (const k of Object.keys(out)) if (k in t) out[k] = t[k];
    for (const k of ["positive","negative","polar","special","hydrophobic","other"]) {
      if (t[k]) {
        out[k].fill = t[k].fill || out[k].fill;
        out[k].stroke = t[k].stroke || out[k].stroke;
      }
    }
    return out;
  } catch (e) {
    return JSON.parse(JSON.stringify(THEME_DEFAULT));
  }
}
let THEME = loadTheme();

const POS_CHARGE = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"], "HIS": ["ND1", "NE2"]};
const NEG_CHARGE = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]};

function isPosCharge(res3, atomName) {
  res3 = (res3||"").toUpperCase();
  atomName = (atomName||"").toUpperCase();
  return (POS_CHARGE[res3] || []).includes(atomName);
}
function isNegCharge(res3, atomName) {
  res3 = (res3||"").toUpperCase();
  atomName = (atomName||"").toUpperCase();
  return (NEG_CHARGE[res3] || []).includes(atomName);
}

function hexToRgb(hex) {
  const h = (hex||"").replace("#","");
  if (h.length !== 6) return [0,0,0];
  return [parseInt(h.slice(0,2),16), parseInt(h.slice(2,4),16), parseInt(h.slice(4,6),16)];
}
function rgba(hex, a) {
  const [r,g,b] = hexToRgb(hex);
  return `rgba(${r},${g},${b},${a})`;
}

function syncColorInputs() {
  if (!themePanel) return;
  $("cPosFill").value = THEME.positive.fill;
  $("cPosStroke").value = THEME.positive.stroke;
  $("cNegFill").value = THEME.negative.fill;
  $("cNegStroke").value = THEME.negative.stroke;
  $("cPolFill").value = THEME.polar.fill;
  $("cPolStroke").value = THEME.polar.stroke;
  $("cSpeFill").value = THEME.special.fill;
  $("cSpeStroke").value = THEME.special.stroke;
  $("cHydFill").value = THEME.hydrophobic.fill;
  $("cHydStroke").value = THEME.hydrophobic.stroke;
  $("cPosAtom").value = THEME.posAtom;
  $("cNegAtom").value = THEME.negAtom;
  $("cDisAtom").value = THEME.disAtom;
}

function applyTheme() {
  for (const [rid, r] of Object.entries(RES)) {
    const g = document.getElementById("res_" + slug(rid));
    if (!g) continue;
    const cls = r.cls || "other";
    const c = THEME[cls] || THEME.other;

    if (r.focus) {
      const center = g.querySelector("circle.center");
      if (center) {
        center.setAttribute("fill", c.fill);
        center.setAttribute("opacity", String(THEME.alpha));
      }
      const t = g.querySelector("text.resLabel");
      if (t) t.setAttribute("fill", c.stroke);
    } else {
      const node = g.querySelector("circle.node");
      const t = g.querySelector("text.nodeLabel");
      if (node) {
        node.setAttribute("fill", c.fill);
        node.setAttribute("fill-opacity", String(THEME.alpha));
        node.setAttribute("stroke", c.stroke);
      }
      if (t) t.setAttribute("fill", c.stroke);
    }
  }

  document.querySelectorAll(".lgPill").forEach(el => {
    const cls = el.getAttribute("data-cls") || "other";
    const c = THEME[cls] || THEME.other;
    el.style.background = rgba(c.fill, THEME.alpha);
    el.style.borderColor = c.stroke;
    el.style.color = c.stroke;
  });

  document.querySelectorAll(".lgDot").forEach(el => {
    const key = el.getAttribute("data-atomkey") || "";
    if (key === "posAtom") el.style.background = THEME.posAtom;
    else if (key === "negAtom") el.style.background = THEME.negAtom;
    else if (key === "disAtom") el.style.background = THEME.disAtom;
  });

  document.querySelectorAll("circle.atom").forEach(el => {
    const name = el.getAttribute("data-name") || "";
    const g = el.closest("g.res");
    const res3 = g ? (g.getAttribute("data-res3") || "") : "";
    if ((res3 || "").toUpperCase() === "CYS" && (name || "").toUpperCase() === "SG") {
      el.setAttribute("fill", THEME.disAtom);
    } else if (isPosCharge(res3, name)) el.setAttribute("fill", THEME.posAtom);
    else if (isNegCharge(res3, name)) el.setAttribute("fill", THEME.negAtom);
  });

  syncColorInputs();
}

function saveTheme() {
  localStorage.setItem("dimplot_theme_v1", JSON.stringify(THEME));
}

function wireColorInputs() {
  const pairs = [
    ["cPosFill", ["positive","fill"]], ["cPosStroke", ["positive","stroke"]],
    ["cNegFill", ["negative","fill"]], ["cNegStroke", ["negative","stroke"]],
    ["cPolFill", ["polar","fill"]], ["cPolStroke", ["polar","stroke"]],
    ["cSpeFill", ["special","fill"]], ["cSpeStroke", ["special","stroke"]],
    ["cHydFill", ["hydrophobic","fill"]], ["cHydStroke", ["hydrophobic","stroke"]],
  ];
  for (const [id, path] of pairs) {
    $(id).addEventListener("input", e => {
      THEME[path[0]][path[1]] = e.target.value;
      saveTheme();
      applyTheme();
    });
  }
  $("cPosAtom").addEventListener("input", e => {
    THEME.posAtom = e.target.value;
    saveTheme();
    applyTheme();
  });
  $("cNegAtom").addEventListener("input", e => {
    THEME.negAtom = e.target.value;
    saveTheme();
    applyTheme();
  });
  $("cDisAtom").addEventListener("input", e => {
    THEME.disAtom = e.target.value;
    saveTheme();
    applyTheme();
  });
}


// -----------------------------
// Chain label renaming (GUI)
// -----------------------------
function applyChainNames(nameA, nameB) {{
  const a = document.getElementById("chainA_label");
  const b = document.getElementById("chainB_label");
  if (a) a.textContent = nameA || "Chain A";
  if (b) b.textContent = nameB || "Chain B";
}}

function loadChainNames() {{
  const nameA = localStorage.getItem("chainAName") || "Chain A";
  const nameB = localStorage.getItem("chainBName") || "Chain B";
  return [nameA, nameB];
}}

function saveChainNames(nameA, nameB) {{
  localStorage.setItem("chainAName", nameA || "Chain A");
  localStorage.setItem("chainBName", nameB || "Chain B");
}}

function wireChainInputs() {{
  const inpA = document.getElementById("chainAName");
  const inpB = document.getElementById("chainBName");
  const btnReset = document.getElementById("btnChainReset");
  if (!inpA || !inpB) return;

  const [nameA0, nameB0] = loadChainNames();
  inpA.value = nameA0;
  inpB.value = nameB0;
  applyChainNames(nameA0, nameB0);

  function commit() {{
    const nameA = (inpA.value || "").trim() || "Chain A";
    const nameB = (inpB.value || "").trim() || "Chain B";
    saveChainNames(nameA, nameB);
    applyChainNames(nameA, nameB);
  }}

  inpA.addEventListener("input", commit);
  inpB.addEventListener("input", commit);

  if (btnReset) {{
    btnReset.addEventListener("click", () => {{
      inpA.value = "Chain A";
      inpB.value = "Chain B";
      commit();
    }});
  }}
}}

let aaMode = (localStorage.getItem("dimplot_aa_mode") || "3");
function applyAAMode() {
  for (const [rid, r] of Object.entries(RES)) {
    const g = document.getElementById("res_" + slug(rid));
    if (!g) continue;
    const txt = r.focus ? g.querySelector("text.resLabel") : g.querySelector("text.nodeLabel");
    if (!txt) continue;
    txt.textContent = (aaMode === "1") ? (r.label1 || r.label3 || r.label) : (r.label3 || r.label);
  }
  btnAA.textContent = (aaMode === "1") ? "AA labels: 1-letter" : "AA labels: 3-letter";
}
function saveAAMode() { localStorage.setItem("dimplot_aa_mode", aaMode); }

btnAA.addEventListener("click", () => {
  aaMode = (aaMode === "1") ? "3" : "1";
  saveAAMode();
  applyAAMode();
});

btnTheme.addEventListener("click", () => {
  const open = (themePanel.style.display === "none");
  themePanel.style.display = open ? "" : "none";
  if (open) {
    positionThemePanel();
    syncColorInputs();
  }
});

window.addEventListener("resize", () => {
  if (themePanel.style.display !== "none") positionThemePanel();
});
window.addEventListener("orientationchange", () => {
  if (themePanel.style.display !== "none") {
    // Wait one frame so viewport sizes settle
    requestAnimationFrame(() => positionThemePanel());
  }
});

btnThemeClose.addEventListener("click", () => {
  themePanel.style.display = "none";
});
btnThemeReset.addEventListener("click", () => {
  THEME = JSON.parse(JSON.stringify(THEME_DEFAULT));
  saveTheme();
  applyTheme();
});

wireColorInputs();

function positionThemePanel() {
  // Keep theme panel within viewport (esp. mobile landscape) and above legend bar.
  try {
    const topbarRect = topbar.getBoundingClientRect();
    const legendRect = legendbar.getBoundingClientRect();
    const vw = window.innerWidth || document.documentElement.clientWidth;
    const vh = window.innerHeight || document.documentElement.clientHeight;

    // Use fixed positioning so the panel never goes off-screen when the page is zoomed/scaled.
    themePanel.style.position = "fixed";
    themePanel.style.right = "10px";

    const top = Math.max(10, Math.round(topbarRect.bottom + 8));
    const safeBottom = Math.min(vh - 10, Math.round(legendRect.top - 8));
    const maxH = Math.max(120, safeBottom - top);

    themePanel.style.top = top + "px";
    themePanel.style.maxHeight = maxH + "px";

    // If extremely constrained, allow the panel to take full height (still scrollable).
    if (maxH < 180) {
      themePanel.style.top = "10px";
      themePanel.style.maxHeight = Math.max(120, vh - 20) + "px";
    }

    // Keep width sensible on narrow screens
    themePanel.style.width = "min(92vw, 380px)";
    themePanel.style.zIndex = "30";
  } catch (e) {}
}

wireChainInputs();


// -----------------------------
// Session export/import (positions + view + options) (added)
// -----------------------------
const btnSessionExport = document.getElementById("btnSessionExport");
const btnSessionImport = document.getElementById("btnSessionImport");
const sessionModal = document.getElementById("sessionModal");
const btnSessionClose = document.getElementById("btnSessionClose");
const btnSessionCopy = document.getElementById("btnSessionCopy");
const btnSessionDownload = document.getElementById("btnSessionDownload");
const btnSessionApply = document.getElementById("btnSessionApply");
const sessionTextarea = document.getElementById("sessionTextarea");
const sessionFile = document.getElementById("sessionFile");

function openSessionModal(text="") {
  if (!sessionModal) return;
  sessionTextarea.value = text || "";
  sessionModal.style.display = "block";
  sessionTextarea.focus();
  sessionTextarea.select();
}

function closeSessionModal() {
  if (!sessionModal) return;
  sessionModal.style.display = "none";
}

function _safeCopy(text) {
  // Try clipboard API first (may fail on file://)
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).catch(() => {});
    return;
  }
  // Fallback
  try {
    sessionTextarea.focus();
    sessionTextarea.select();
    document.execCommand("copy");
  } catch (e) {}
}

function collectSessionState() {
  // capture residue positions + focus rotations/atom offsets
  const residues = {};
  for (const [rid, r] of Object.entries(RES)) {
    residues[rid] = {
      x: r.x, y: r.y,
      rot: (r.rot || 0),
      aox: (r.aox || 0), aoy: (r.aoy || 0)
    };
  }
  // chain labels
  let chainA = (localStorage.getItem("dimplot_chainA") || "Chain A");
  let chainB = (localStorage.getItem("dimplot_chainB") || "Chain B");
  try {
    const a = document.getElementById("chainAName");
    const b = document.getElementById("chainBName");
    if (a && a.value) chainA = a.value.trim() || chainA;
    if (b && b.value) chainB = b.value.trim() || chainB;
  } catch(e) {}

  return {
    schema: "DIMPLOTViewer.session.v1",
    createdAt: new Date().toISOString(),
    theme: THEME,
    chainA, chainB,
    aaMode: (typeof aaMode !== "undefined" ? aaMode : "3"),
    view: {x: view.x, y: view.y, k: view.k},
    flags: {
      dragEnabled, chainSideLock, panZoomEnabled,
      contactsVisible, hbOnly
    },
    residues
  };
}

function applySessionState(st) {
  if (!st || typeof st !== "object") return;

  // Theme
  if (st.theme && typeof st.theme === "object") {
    for (const k of Object.keys(THEME)) {
      if (k in st.theme) THEME[k] = st.theme[k];
    }
    applyTheme();
    // update color inputs if present
    try { wireColorInputs();

function positionThemePanel() {
  // Keep theme panel within viewport (esp. mobile landscape) and above legend bar.
  try {
    const topbarRect = topbar.getBoundingClientRect();
    const legendRect = legendbar.getBoundingClientRect();
    const vw = window.innerWidth || document.documentElement.clientWidth;
    const vh = window.innerHeight || document.documentElement.clientHeight;

    // Use fixed positioning so the panel never goes off-screen when the page is zoomed/scaled.
    themePanel.style.position = "fixed";
    themePanel.style.right = "10px";

    const top = Math.max(10, Math.round(topbarRect.bottom + 8));
    const safeBottom = Math.min(vh - 10, Math.round(legendRect.top - 8));
    const maxH = Math.max(120, safeBottom - top);

    themePanel.style.top = top + "px";
    themePanel.style.maxHeight = maxH + "px";

    // If extremely constrained, allow the panel to take full height (still scrollable).
    if (maxH < 180) {
      themePanel.style.top = "10px";
      themePanel.style.maxHeight = Math.max(120, vh - 20) + "px";
    }

    // Keep width sensible on narrow screens
    themePanel.style.width = "min(92vw, 380px)";
    themePanel.style.zIndex = "30";
  } catch (e) {}
}
 } catch(e) {}
  }

  // Chain labels
  if (typeof st.chainA === "string" || typeof st.chainB === "string") {
    const a = (st.chainA || "Chain A");
    const b = (st.chainB || "Chain B");
    try {
      localStorage.setItem("dimplot_chainA", a);
      localStorage.setItem("dimplot_chainB", b);
    } catch(e) {}
    try {
      const inpA = document.getElementById("chainAName");
      const inpB = document.getElementById("chainBName");
      if (inpA) inpA.value = a;
      if (inpB) inpB.value = b;
      applyChainNames(a,b);
    } catch(e) {}
  }

  // AA mode
  if (st.aaMode === "1" || st.aaMode === "3") {
    aaMode = st.aaMode;
    try { localStorage.setItem("dimplot_aa_mode", aaMode); } catch(e) {}
    applyAAMode();
    try {
      const btn = document.getElementById("btnAA");
      if (btn) btn.textContent = (aaMode === "1" ? "1-letter" : "3-letter");
    } catch(e) {}
  }

  // View (zoom/pan)
  if (st.view && typeof st.view === "object") {
    if (typeof st.view.x === "number") view.x = st.view.x;
    if (typeof st.view.y === "number") view.y = st.view.y;
    if (typeof st.view.k === "number") view.k = st.view.k;
    applyView();
  }

  // Flags
  if (st.flags && typeof st.flags === "object") {
    if (typeof st.flags.dragEnabled === "boolean") dragEnabled = st.flags.dragEnabled;
    if (typeof st.flags.chainSideLock === "boolean") chainSideLock = st.flags.chainSideLock;
    if (typeof st.flags.panZoomEnabled === "boolean") panZoomEnabled = st.flags.panZoomEnabled;
    if (typeof st.flags.contactsVisible === "boolean") contactsVisible = st.flags.contactsVisible;
    if (typeof st.flags.hbOnly === "boolean") hbOnly = st.flags.hbOnly;

    // sync buttons (do NOT trigger side effects beyond display)
    try {
      btnDrag.classList.toggle('on', dragEnabled);
      btnDrag.textContent = dragEnabled ? "Drag ON" : "Drag OFF";
    } catch(e) {}
    try {
      btnLock.classList.toggle('on', chainSideLock);
      btnLock.textContent = chainSideLock ? "Lock ON" : "Lock OFF";
    } catch(e) {}
    try {
      btnPan.classList.toggle('on', panZoomEnabled);
      btnPan.textContent = panZoomEnabled ? "Pan ON" : "Pan OFF";
    } catch(e) {}
    try {
      btnContacts.classList.toggle('on', contactsVisible);
      btnContacts.textContent = contactsVisible ? "Contacts ON" : "Contacts OFF";
    } catch(e) {}
    try {
      btnHBOnly.classList.toggle('on', hbOnly);
      btnHBOnly.textContent = hbOnly ? "H-bond only: ON" : "H-bond only: OFF";
    } catch(e) {}
  }

  // Residues positions/rotations
  if (st.residues && typeof st.residues === "object") {
    for (const [rid, s] of Object.entries(st.residues)) {
      const r = RES[rid];
      if (!r || !s) continue;
      if (typeof s.x === "number") r.x = s.x;
      if (typeof s.y === "number") r.y = s.y;
      if (typeof s.rot === "number") r.rot = s.rot;
      if (typeof s.aox === "number") r.aox = s.aox;
      if (typeof s.aoy === "number") r.aoy = s.aoy;

      const g = document.getElementById("res_" + slug(rid));
      if (g) g.setAttribute("transform", `translate(${r.x},${r.y})`);
      if (r.focus) setAtomBlockRotation(rid);
    }
    updateAllEdges();
  }

  // Re-apply filters/visibility and HB top labels
  try { applyResidueFilter(); } catch(e) {}
  try { applyHbondVisibilityByFilter(); } catch(e) {}
  try { autoPlaceFocusTopLabelsAndPull(); } catch(e) {}
  try { updateAllEdges(); } catch(e) {}
}

function downloadText(filename, text) {
  try {
    const blob = new Blob([text], {type:"application/json"});
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1500);
  } catch (e) {}
}

if (btnSessionExport) {
  btnSessionExport.addEventListener("click", (e) => {
    e.preventDefault();
    e.stopPropagation();
    const st = collectSessionState();
    const json = JSON.stringify(st, null, 2);
    openSessionModal(json);
  });
}
if (btnSessionImport) {
  btnSessionImport.addEventListener("click", (e) => {
    e.preventDefault();
    e.stopPropagation();
    openSessionModal("");
    // open file picker as convenience (user can cancel)
    try { sessionFile.click(); } catch(err) {}
  });
}
if (btnSessionClose) btnSessionClose.addEventListener("click", (e)=>{ e.preventDefault(); closeSessionModal(); });
if (sessionModal) {
  // clicking backdrop closes
  sessionModal.addEventListener("pointerdown", (e)=>{ if (e.target === sessionModal) closeSessionModal(); });
  // prevent canvas gesture stealing
  sessionModal.addEventListener("pointerdown", (e)=>{ e.stopPropagation(); }, true);
}

if (btnSessionCopy) btnSessionCopy.addEventListener("click", ()=>{ _safeCopy(sessionTextarea.value || ""); });

if (btnSessionDownload) btnSessionDownload.addEventListener("click", ()=>{
  const txt = sessionTextarea.value || "";
  if (!txt.trim()) return;
  downloadText("2DPPIViewer_session.json", txt);
});

if (btnSessionApply) btnSessionApply.addEventListener("click", ()=>{
  const txt = sessionTextarea.value || "";
  if (!txt.trim()) return;
  try {
    const st = JSON.parse(txt);
    applySessionState(st);
    closeSessionModal();
  } catch (e) {
    alert("Invalid JSON.");
  }
});

if (sessionFile) {
  sessionFile.addEventListener("change", async ()=>{
    const f = sessionFile.files && sessionFile.files[0];
    if (!f) return;
    try {
      const txt = await f.text();
      sessionTextarea.value = txt;
    } catch(e) {}
  });
}

// Prevent theme panel inputs from triggering canvas drag on mobile
if (themePanel) {
  themePanel.addEventListener("pointerdown", (e)=>{ e.stopPropagation(); }, true);
  themePanel.addEventListener("touchstart", (e)=>{ e.stopPropagation(); }, true);
}


// -----------------------------
// Touch pinch-zoom (mobile) (added)
// - Only intercept when 2 fingers are active, so single-finger drag stays unchanged.
// -----------------------------
(function(){
  const pts = new Map(); // pointerId -> {x,y}
  let pinching = false;
  let pinch0 = null; // {dist, midx, midy, viewx, viewy, viewk}

  function dist(a,b){ const dx=a.x-b.x, dy=a.y-b.y; return Math.hypot(dx,dy); }
  function mid(a,b){ return {x:(a.x+b.x)/2, y:(a.y+b.y)/2}; }

  svg.addEventListener("pointerdown", (evt)=>{
    if (evt.pointerType !== "touch") return;
    if (evt.target.closest("#topbar") || evt.target.closest("#themePanel") || evt.target.closest("#sessionModal")) return;
    const p = svgPoint(evt);
    pts.set(evt.pointerId, {x:p.x, y:p.y});
    if (pts.size === 2) {
      const arr = Array.from(pts.values());
      const d0 = dist(arr[0], arr[1]);
      const m0 = mid(arr[0], arr[1]);
      pinching = true;
      pinch0 = {dist:d0, midx:m0.x, midy:m0.y, viewx:view.x, viewy:view.y, viewk:view.k};
      // Stop other handlers (drag/pan) only once pinch starts
      evt.preventDefault();
      evt.stopImmediatePropagation();
    }
  }, {capture:true});

  svg.addEventListener("pointermove", (evt)=>{
    if (!pinching) return;
    if (evt.pointerType !== "touch") return;
    if (!pts.has(evt.pointerId)) return;
    const p = svgPoint(evt);
    pts.set(evt.pointerId, {x:p.x, y:p.y});
    if (pts.size !== 2 || !pinch0) return;

    const arr = Array.from(pts.values());
    const d = dist(arr[0], arr[1]);
    const m = mid(arr[0], arr[1]);
    const scale = (d && pinch0.dist) ? (d / pinch0.dist) : 1;

    const newK = Math.min(8, Math.max(0.2, pinch0.viewk * scale));
    view.k = newK;
    // Pan with midpoint movement (screen space)
    view.x = pinch0.viewx + (m.x - pinch0.midx);
    view.y = pinch0.viewy + (m.y - pinch0.midy);

    applyView();
    evt.preventDefault();
    evt.stopImmediatePropagation();
  }, {capture:true});

  function end(evt){
    if (evt && pts.has(evt.pointerId)) pts.delete(evt.pointerId);
    if (pts.size < 2) { pinching = false; pinch0 = null; }
  }
  svg.addEventListener("pointerup",   (evt)=>{ end(evt); }, {capture:true});
  svg.addEventListener("pointercancel",(evt)=>{ end(evt); }, {capture:true});
})();
"""

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{escape(UI["title"])}</title>
<style>
  html, body {{
    margin:0; height:100%; overflow:hidden;
    font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif;
  }}

  #topbar {{
    position: fixed; left: 10px; top: 10px; z-index: 10;
    display:flex; gap:10px; align-items:center; flex-wrap:wrap;
    max-width: calc(100% - 20px);
  }}

  .btn {{
    border:1px solid #c7c7c7; background:#fff; padding:8px 12px; border-radius:10px;
    font-weight:800; cursor:pointer; user-select:none;
    box-shadow: 0 1px 6px rgba(0,0,0,0.08);
  }}
  .btn.on {{ background:#e8f0fe; border-color:#7aa5ff; }}

  #hint {{
    margin-left: 10px; color:#444; font-size: 13px;
    max-width: 980px;
    line-height: 1.25;
  }}

  #legendbar {{
    position: fixed;
    left: 10px;
    right: 10px;
    bottom: 10px;
    z-index: 10;
    display:flex;
    gap: 12px;
    align-items:center;
    justify-content: center;
    flex-wrap: wrap;
    background: rgba(255,255,255,0.92);
    border: 1px solid #d6d6d6;
    border-radius: 12px;
    padding: 8px 10px;
    box-shadow: 0 1px 6px rgba(0,0,0,0.08);
    font-size: 12.5px;
    color: #222;
  }}
  .legendGroup {{ display:flex; align-items:center; gap:10px; flex-wrap:wrap; }}
  .legendSpacer {{ width: 18px; height: 18px; border-left: 1px solid rgba(0,0,0,0.18); margin: 0 6px; }}

  /* Residue legend pills */
  .lgPill {{
    display:inline-flex; align-items:center; justify-content:center;
    padding: 6px 10px;
    border-radius: 999px;
    border: 2px solid rgba(0,0,0,0.25);
    font-weight: 800;
    letter-spacing: 0.1px;
    user-select:none;
  }}

  /* Atom legend (dots) */
  .lgAtomItem {{ display:inline-flex; align-items:center; gap:8px; font-weight:800; }}
  .lgDot {{
    width: 12px; height: 12px; border-radius: 999px;
    border: 1px solid rgba(0,0,0,0.28);
    display:inline-block;
  }}
  .lgAtomText {{ color: #222; }}

  svg {{ width:100%; height:100%; background:#ffffff; overflow: visible; }}

  .edgeLabel, .nodeLabel, .resLabel, .atomLabel {{
    pointer-events:none;
    stroke: none;
  }}

  .res {{ cursor: grab; }}
  .res:active {{ cursor: grabbing; }}

  /* Mobile + touch improvements */
  #svg {{ touch-action: none; }}
  #themePanel input, #themePanel textarea, #themePanel select {{ font-size:16px; }}
  #themePanel {{ -webkit-overflow-scrolling: touch; }}

  @media (max-width: 520px) {{
    #topbar {{ left: 6px; top: 6px; gap:6px; }}
    .btn {{ padding:6px 9px; border-radius:10px; font-size:12px; }}
    #themePanel {{ width: min(92vw, 340px); max-height: calc(100vh - 120px); }}
    /* Make SVG labels smaller on phones */
    text.resLabel {{ font-size: 10px; }}
    text.atomLabel {{ font-size: 9px; }}
    text.chainLabel {{ font-size: 12px; }}
    #legendbar {{ font-size: 12px; }}
  }}


  /* Theme panel sizing + stacking */
  #themePanel {{ overflow:auto; max-height: calc(100vh - 140px); z-index: 30; }}

  /* Landscape phone: reduce UI footprint */
  @media (orientation: landscape) and (max-height: 520px) {{
    #hint {{ display:none; }}
    #topbar {{ gap:6px; }}
    .btn {{ padding:6px 9px; border-radius:10px; font-size:12px; }}
    #legendbar {{ padding: 4px 6px; gap: 8px; font-size: 11.5px; max-height: 72px; overflow:auto; }}
    .lgPill {{ padding: 4px 8px; }}
    .lgDot {{ width: 10px; height: 10px; }}
    #themePanel {{ width: min(70vw, 380px); max-height: calc(100vh - 90px); }}
    text.resLabel {{ font-size: 9px; }}
    text.atomLabel {{ font-size: 8px; }}
    text.chainLabel {{ font-size: 11px; }}
  }}

</style>
</head>
<body>

<div id="topbar">
  <div class="btn on" id="btnDrag">{escape(UI["btn_drag_on"])}</div>
  <div class="btn on" id="btnLock">{escape(UI["btn_lock_on"])}</div>
  <div class="btn on" id="btnPan">{escape(UI["btn_pan_on"])}</div>
  <div class="btn on" id="btnContacts">{escape(UI["btn_contacts_on"])}</div>
  <div class="btn" id="btnHBOnly">{escape(UI["btn_hbonly_off"])}</div>
  <div class="btn" id="btnAA">AA labels: 3-letter</div>
  
  <div class="btn" id="btnFit">{escape(UI["btn_fit"])}</div>
  <div class="btn" id="btnReset">{escape(UI["btn_reset"])}</div>
  <div class="btn" id="btnExport">{escape(UI["btn_export"])}</div>
  <div id="themeWrap" style="position:relative; display:inline-block;">
  <div class="btn" id="btnTheme">Theme</div>
  <div id="themePanel" style="display:none; position:absolute; top: calc(100% + 8px); right:0; padding:10px 12px; background:rgba(255,255,255,0.96); border:1px solid #d6d6d6; border-radius:12px; box-shadow:0 6px 18px rgba(0,0,0,0.16); z-index:999; width: 320px; max-height: calc(100vh - 140px); overflow:auto;">

    <div style="font-weight:900; margin-bottom:6px;">Colors</div>
    <div style="display:grid; grid-template-columns: 110px 90px 90px; gap:6px 10px; align-items:center; font-size:12px;">
      <div></div><div style="font-weight:800;">Fill</div><div style="font-weight:800;">Stroke</div>
      <div>Positive</div><input type="color" id="cPosFill"><input type="color" id="cPosStroke">
      <div>Negative</div><input type="color" id="cNegFill"><input type="color" id="cNegStroke">
      <div>Polar</div><input type="color" id="cPolFill"><input type="color" id="cPolStroke">
      <div>Special</div><input type="color" id="cSpeFill"><input type="color" id="cSpeStroke">
      <div>Hydrophobic</div><input type="color" id="cHydFill"><input type="color" id="cHydStroke">
      <div style="grid-column:1 / span 3; height:1px; background:rgba(0,0,0,0.15); margin:4px 0;"></div>
      <div>+charge atom</div><input type="color" id="cPosAtom"><div></div>
      <div>-charge atom</div><input type="color" id="cNegAtom"><div></div>
      <div>Disulphide atom</div><input type="color" id="cDisAtom"><div></div>
    </div>

    <div style="height:10px;"></div>
    <div style="font-weight:900; margin:4px 0 6px;">Chain labels</div>
    <div style="display:grid; grid-template-columns: 90px 1fr; gap:6px 10px; align-items:center; font-size:12px;">
      <div>Chain A</div><input type="text" id="chainAName" placeholder="Chain A" style="width:100%; padding:6px 8px; border:1px solid #d6d6d6; border-radius:8px; font-size:12px;">
      <div>Chain B</div><input type="text" id="chainBName" placeholder="Chain B" style="width:100%; padding:6px 8px; border:1px solid #d6d6d6; border-radius:8px; font-size:12px;">
    </div>

    <hr style="border:none; border-top:1px solid #e5e5e5; margin:12px 0 10px;">
    <div style="font-weight:900; margin-bottom:6px;">Session</div>
    <div style="display:flex; gap:8px; flex-wrap:wrap;">
      <div class="btn" id="btnSessionExport" style="padding:6px 10px;">Export JSON</div>
      <div class="btn" id="btnSessionImport" style="padding:6px 10px;">Import JSON</div>
    </div>
    <div style="font-size:12px; opacity:0.75; margin-top:6px; line-height:1.25;">
      Exports/imports: colors, chain names, AA label mode, view (zoom/pan), residue positions/rotations.
    </div>

    <div style="display:flex; gap:8px; margin-top:8px; justify-content:flex-end;">
      <div class="btn" id="btnChainReset" style="padding:6px 10px;">Reset chains</div>
    </div>
    <div style="display:flex; gap:8px; margin-top:8px; justify-content:flex-end;">
      <div class="btn" id="btnThemeReset" style="padding:6px 10px;">Reset colors</div>
      <div class="btn on" id="btnThemeClose" style="padding:6px 10px;">Close</div>
    </div>
  </div>
</div>

  <div id="hint">{escape(UI["hint"])}</div>
</div>


<div id="sessionModal" style="display:none; position:fixed; inset:0; z-index:10000; background:rgba(0,0,0,0.45);">
  <div style="position:absolute; left:50%; top:50%; transform:translate(-50%,-50%); width:min(92vw, 720px); background:rgba(255,255,255,0.98); border:1px solid #d6d6d6; border-radius:14px; box-shadow:0 10px 28px rgba(0,0,0,0.25); padding:12px;">
    <div style="display:flex; justify-content:space-between; align-items:center; gap:10px; margin-bottom:8px;">
      <div style="font-weight:900;">Session JSON</div>
      <div class="btn" id="btnSessionClose" style="padding:6px 10px;">Close</div>
    </div>
    <textarea id="sessionTextarea" spellcheck="false" style="width:100%; height: min(55vh, 420px); resize:vertical; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace; font-size:12px; padding:10px; border:1px solid #ddd; border-radius:10px;"></textarea>
    <div style="display:flex; gap:8px; flex-wrap:wrap; justify-content:flex-end; margin-top:8px;">
      <input id="sessionFile" type="file" accept="application/json" style="display:none;">
      <div class="btn" id="btnSessionCopy" style="padding:6px 10px;">Copy</div>
      <div class="btn" id="btnSessionDownload" style="padding:6px 10px;">Download</div>
      <div class="btn on" id="btnSessionApply" style="padding:6px 10px;">Apply</div>
    </div>
    <div style="font-size:12px; opacity:0.75; margin-top:6px; line-height:1.25;">
      Tip: On mobile, use Copy and paste into Notes/email to move to PC, or Download if supported.
    </div>
  </div>
</div>


<div id="legendbar">
  {legend_html}
</div>

<svg id="svg" viewBox="0 0 {width} {height}">
  <g id="viewport">
  {svg_str}
  </g>
</svg>

<script>
const W = {width}, H = {height};
const BOUNDARY_Y = {boundary_y:.3f};

const RES = {json.dumps(res_data, ensure_ascii=False)};
// per-focus atom block offset (to pull toward boundary when label would overlap another residue)
for (const [rid, r] of Object.entries(RES)) {{
  if (r.aox === undefined) r.aox = 0;
  if (r.aoy === undefined) r.aoy = 0;
  r._aox0 = r.aox;
  r._aoy0 = r.aoy;
}}
const LABEL_TOP_GAP = 10;   // gap between atom diagram and residue label (larger => label higher)
const PULL_STEP = 4;        // px per iteration when pulling toward boundary
const PULL_MAX_ITERS = 30;  // safety

const EDGES = {json.dumps(edges, ensure_ascii=False)};

const SIMPLE_R = {simple_r:.3f};
const FOCUS_R  = {focus_r:.3f};

const RELAX_ITERS = {relax_iters};
const RELAX_GAP   = {relax_gap_px};
const RELAX_STEP  = {relax_step};

const RELAX_ON_POINTERUP = {str(bool(relax_on_pointerup)).lower()};

let dragEnabled = true;
let chainSideLock = true;
let panZoomEnabled = true;
let contactsVisible = true;

// New: hide residues that are not in H-bond network
let hbOnly = false;

const btnDrag = document.getElementById('btnDrag');
const btnLock = document.getElementById('btnLock');
const btnPan  = document.getElementById('btnPan');
const btnContacts = document.getElementById('btnContacts');
const btnHBOnly = document.getElementById('btnHBOnly');
const btnFit  = document.getElementById('btnFit');
const btnReset= document.getElementById('btnReset');
const btnExport= document.getElementById('btnExport');

const svg = document.getElementById('svg');
const viewport = document.getElementById('viewport');

function slug(s) {{ return s.replace(/[^A-Za-z0-9_]+/g, '_'); }}

for (const [rid, r] of Object.entries(RES)) {{
  r._x0 = r.x; r._y0 = r.y;
  r.rot = 0;
  r._rot0 = 0;
}}

// Build set of residues participating in ANY H-bond
const HB_RES = new Set();
for (const e of EDGES) {{
  if (e.type === 'hbond') {{
    HB_RES.add(e.ra);
    HB_RES.add(e.rb);
  }}
}}

btnDrag.addEventListener('click', () => {{
  dragEnabled = !dragEnabled;
  btnDrag.classList.toggle('on', dragEnabled);
  btnDrag.textContent = dragEnabled ? {json.dumps(UI["btn_drag_on"])} : {json.dumps(UI["btn_drag_off"])};
}});

btnLock.addEventListener('click', () => {{
  chainSideLock = !chainSideLock;
  btnLock.classList.toggle('on', chainSideLock);
  btnLock.textContent = chainSideLock ? {json.dumps(UI["btn_lock_on"])} : {json.dumps(UI["btn_lock_off"])};
  resolveOverlaps();
  updateAllEdges();
  autoPlaceFocusTopLabelsAndPull();
}});

btnPan.addEventListener('click', () => {{
  panZoomEnabled = !panZoomEnabled;
  btnPan.classList.toggle('on', panZoomEnabled);
  btnPan.textContent = panZoomEnabled ? {json.dumps(UI["btn_pan_on"])} : {json.dumps(UI["btn_pan_off"])};
}});

function applyContactsVisibility() {{
  const els = document.querySelectorAll('line.edge.contact');
  els.forEach(el => {{
    el.style.display = contactsVisible ? '' : 'none';
  }});
}}

function applyResidueFilter() {{
  // hbOnly: show only residues that are in the H-bond network
  for (const rid of Object.keys(RES)) {{
    const g = document.getElementById('res_' + slug(rid));
    if (!g) continue;
    if (!hbOnly) {{
      g.style.display = '';
    }} else {{
      g.style.display = HB_RES.has(rid) ? '' : 'none';
    }}
  }}
  applyContactsVisibility();

  // Also hide H-bond / disulphide lines if either endpoint residue is hidden (safety)
  for (const e of EDGES) {{
    if (e.type !== 'hbond' && e.type !== 'disulphide') continue;
    const prefix = (e.type === 'hbond') ? 'h_' : 'd_';
    const lid = prefix + slug(e.ra) + '_' + e.a + '__' + slug(e.rb) + '_' + e.b;
    const line = document.getElementById(lid);
    const text = document.getElementById(lid + '_t');

    const show = (!hbOnly) || (HB_RES.has(e.ra) && HB_RES.has(e.rb));
    if (line) line.style.display = show ? '' : 'none';
    if (text) text.style.display = show ? '' : 'none';
  }}

  updateAllEdges();
}}

btnContacts.addEventListener('click', () => {{
  const turningOn = !contactsVisible;
  contactsVisible = !contactsVisible;
  btnContacts.classList.toggle('on', contactsVisible);
  btnContacts.textContent = contactsVisible ? {json.dumps(UI["btn_contacts_on"])} : {json.dumps(UI["btn_contacts_off"])};

  // If user turns ON contacts while in HB-only mode, automatically show all residues
  // (prevents "lines-only" awkwardness).
  if (turningOn && hbOnly) {{
    hbOnly = false;
    btnHBOnly.classList.toggle('on', hbOnly);
    btnHBOnly.textContent = hbOnly ? {json.dumps(UI["btn_hbonly_on"])} : {json.dumps(UI["btn_hbonly_off"])};
    applyResidueFilter();
  }}

  applyContactsVisibility();
  updateAllEdges();
}});

btnHBOnly.addEventListener('click', () => {{
  hbOnly = !hbOnly;
  btnHBOnly.classList.toggle('on', hbOnly);
  btnHBOnly.textContent = hbOnly ? {json.dumps(UI["btn_hbonly_on"])} : {json.dumps(UI["btn_hbonly_off"])};
  applyResidueFilter();
}});

btnFit.addEventListener('click', () => {{
  fitViewToResidues();
}});

function setAtomBlockRotation(rid) {{
  const r = RES[rid];
  if (!r || !r.focus) return;
  const ab = document.getElementById('atomblock_' + slug(rid));
  if (!ab) return;
  const deg = (r.rot || 0);
  ab.setAttribute('transform', `rotate(${{deg}})`);
}}


function _rectOverlap(a, b, pad=0) {{
  return !((a.right + pad) < b.left || (b.right + pad) < a.left || (a.bottom + pad) < b.top || (b.bottom + pad) < a.top);
}}

// Place focus residue label above its atom diagram; if it overlaps another residue node circle,
// pull the atom diagram + label toward the boundary line (keeps initial orientation stable).
function autoPlaceFocusTopLabelsAndPull() {{
  // Build list of other residue circles (simple nodes + focus centers)
  const circles = [];
  for (const [rid2, r2] of Object.entries(RES)) {{
    const g2 = document.getElementById('res_' + slug(rid2));
    if (!g2 || g2.style.display === 'none') continue;
    if (hbOnly && !HB_RES.has(rid2)) continue;

    let c = null;
    if (r2.focus) c = g2.querySelector('circle.center');
    else c = g2.querySelector('circle.node');
    if (!c) continue;
    circles.push({{rid: rid2, el: c}});
  }}

  for (const [rid, r] of Object.entries(RES)) {{
    if (!r || !r.focus) continue;

    const g = document.getElementById('res_' + slug(rid));
    if (!g || g.style.display === 'none') continue;
    if (hbOnly && !HB_RES.has(rid)) continue;

    const label = g.querySelector('text.resLabel');
    const ab = document.getElementById('atomblock_' + slug(rid));
    if (!label || !ab) continue;

    // compute desired label position in local coords (accounting for aox/aoy since bbox ignores transforms)
    let bb;
    try {{ bb = ab.getBBox(); }} catch(e) {{ continue; }}

    const lx = (bb.x + bb.width/2) + (r.aox||0);
    const ly = (bb.y) + (r.aoy||0) - LABEL_TOP_GAP;

    label.setAttribute('text-anchor', 'middle');
    label.setAttribute('dominant-baseline', 'alphabetic');
    label.setAttribute('x', lx);
    label.setAttribute('y', ly);

    // If the label overlaps another residue circle, pull toward boundary by shifting aoy.
    let it = 0;
    while (it < PULL_MAX_ITERS) {{
      const lr = label.getBoundingClientRect();
      let hit = false;

      for (const c of circles) {{
        if (c.rid === rid) continue;
        const cr = c.el.getBoundingClientRect();
        if (_rectOverlap(lr, cr, 2)) {{ hit = true; break; }}
      }}

      if (!hit) break;

      const dir = (r.chain === 'A') ? +1 : (r.chain === 'B' ? -1 : 0);
      if (dir === 0) break;

      r.aoy = (r.aoy||0) + dir * PULL_STEP;

      // apply updated transform and recompute label position
      setAtomBlockRotation(rid);
      const ly2 = (bb.y) + (r.aoy||0) - LABEL_TOP_GAP;
      label.setAttribute('y', ly2);

      it++;
    }}
  }}
}}


btnReset.addEventListener('click', () => {{
  for (const [rid, r] of Object.entries(RES)) {{
    r.x = r._x0; r.y = r._y0;
    r.rot = r._rot0 || 0;
    r.aox = r._aox0 || 0;
    r.aoy = r._aoy0 || 0;
    const g = document.getElementById('res_' + slug(rid));
    if (g) g.setAttribute('transform', `translate(${{r.x}},${{r.y}})`);
    setAtomBlockRotation(rid);
  }}

  view = {{x:0, y:0, k:1}};
  applyView();

  // restore toggles
  dragEnabled = true;
  chainSideLock = true;
  panZoomEnabled = true;
  contactsVisible = true;
  hbOnly = false;

  btnDrag.classList.toggle('on', true);
  btnDrag.textContent = {json.dumps(UI["btn_drag_on"])};

  btnLock.classList.toggle('on', true);
  btnLock.textContent = {json.dumps(UI["btn_lock_on"])};

  btnPan.classList.toggle('on', true);
  btnPan.textContent = {json.dumps(UI["btn_pan_on"])};

  btnContacts.classList.toggle('on', true);
  btnContacts.textContent = {json.dumps(UI["btn_contacts_on"])};

  btnHBOnly.classList.toggle('on', false);
  btnHBOnly.textContent = {json.dumps(UI["btn_hbonly_off"])};

  resolveOverlaps();
  updateAllEdges();
  autoPlaceFocusTopLabelsAndPull();
  fitViewToResidues();
  applyContactsVisibility();
  applyResidueFilter();
}});

function exportCurrentSVG() {{
  const src = document.getElementById('svg');
  const clone = src.cloneNode(true);

  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');

  const PAD = 240;
  clone.setAttribute('viewBox', (-PAD) + ' ' + (-PAD) + ' ' + (W + 2*PAD) + ' ' + (H + 2*PAD));
  clone.setAttribute('width', String(W + 2*PAD));
  clone.setAttribute('height', String(H + 2*PAD));

  const serializer = new XMLSerializer();
  let svgText = serializer.serializeToString(clone);
  if (!svgText.startsWith('<?xml')) {{
    svgText = '<?xml version="1.0" encoding="UTF-8"?>\\n' + svgText;
  }}

  const blob = new Blob([svgText], {{type: 'image/svg+xml;charset=utf-8'}});
  const url = URL.createObjectURL(blob);

  const a = document.createElement('a');
  a.href = url;
  a.download = '2DPPIViewer.svg';
  document.body.appendChild(a);
  a.click();
  a.remove();

  setTimeout(() => URL.revokeObjectURL(url), 500);
}}

btnExport.addEventListener('click', () => {{
  exportCurrentSVG();
}});

function setLine(id, x1,y1,x2,y2) {{
  const el = document.getElementById(id);
  if (!el) return;
  el.setAttribute('x1', x1); el.setAttribute('y1', y1);
  el.setAttribute('x2', x2); el.setAttribute('y2', y2);
}}
function setText(id, x,y, txt) {{
  const el = document.getElementById(id);
  if (!el) return;
  el.setAttribute('x', x); el.setAttribute('y', y);
  el.textContent = txt;
}}

function clampY(chain, y) {{
  if (!chainSideLock) return y;
  const gap = 26;
  if (chain === 'A') return Math.min(y, BOUNDARY_Y - gap);
  if (chain === 'B') return Math.max(y, BOUNDARY_Y + gap);
  return y;
}}

function getAtomAbs(resId, atomIdxStr) {{
  const r = RES[resId];
  if (!r) return null;

  if (r.atoms && r.atoms[atomIdxStr]) {{
    const a = r.atoms[atomIdxStr];
    const ax = a.dx, ay = a.dy;

    const rad = (r.rot || 0) * Math.PI / 180.0;
    const c = Math.cos(rad), s = Math.sin(rad);
    const rx = ax * c - ay * s;
    const ry = ax * s + ay * c;

    return {{ x: r.x + rx, y: r.y + ry }};
  }}

  return {{ x: r.x, y: r.y }};
}}

function updateAllEdges() {{
  // contacts
  for (const e of EDGES) {{
    if (e.type === 'contact') {{
      const u = RES[e.u], v = RES[e.v];
      if (!u || !v) continue;

      const lid = 'c_' + slug(e.u) + '__' + slug(e.v);
      const line = document.getElementById(lid);
      if (!line) continue;

      // In hbOnly mode, show contacts only between visible (HB) residues.
      const show = contactsVisible && (!hbOnly || (HB_RES.has(e.u) && HB_RES.has(e.v)));
      line.style.display = show ? '' : 'none';

      setLine(lid, u.x, u.y, v.x, v.y);


    }}
  }}

  // hbonds / disulphides
  for (const e of EDGES) {{
    if (e.type === 'hbond' || e.type === 'disulphide') {{
      if (hbOnly && (!HB_RES.has(e.ra) || !HB_RES.has(e.rb))) continue;

      const p1 = getAtomAbs(e.ra, e.a);
      const p2 = getAtomAbs(e.rb, e.b);
      if (!p1 || !p2) continue;

      const prefix = (e.type === 'hbond') ? 'h_' : 'd_';
      const lid = prefix + slug(e.ra) + '_' + e.a + '__' + slug(e.rb) + '_' + e.b;
      setLine(lid, p1.x, p1.y, p2.x, p2.y);

      if (e.type === 'hbond') {{
        const mx = (p1.x + p2.x) / 2, my = (p1.y + p2.y) / 2;
        const dx = p2.x - p1.x, dy = p2.y - p1.y;
        const L = Math.hypot(dx, dy) || 1;
        const nx = -dy / L, ny = dx / L;
        const off = 12;
        const tx = mx + nx * off, ty = my + ny * off;

        const distTxt = (e.dist !== null && e.dist !== undefined)
          ? (Number(e.dist).toFixed(2) + ' Å') : '';
        setText(lid + '_t', tx, ty, distTxt);
      }}
    }}
  }}
}}

function nodeRadius(r) {{
  return r.focus ? FOCUS_R : SIMPLE_R;
}}

function resolveOverlaps() {{
  const ids = Object.keys(RES);
  for (let it = 0; it < RELAX_ITERS; it++) {{
    let moved = 0;
    for (let i = 0; i < ids.length; i++) {{
      const A = RES[ids[i]];
      if (!A) continue;

      for (let j = i + 1; j < ids.length; j++) {{
        const B = RES[ids[j]];
        if (!B) continue;

        if (A.chain !== B.chain) continue;

        // if hbOnly: only relax among visible residues
        if (hbOnly && (!HB_RES.has(ids[i]) || !HB_RES.has(ids[j]))) continue;

        const minD = nodeRadius(A) + nodeRadius(B) + RELAX_GAP;
        let dx = B.x - A.x;
        let dy = B.y - A.y;
        let d = Math.hypot(dx, dy);

        if (d < 1e-6) {{ dx = 1; dy = 0; d = 1; }}

        if (d < minD) {{
          const push = (minD - d) * 0.5 * RELAX_STEP;
          const ux = dx / d;
          const uy = dy / d;

          A.x -= ux * push; A.y -= uy * push;
          B.x += ux * push; B.y += uy * push;

          A.y = clampY(A.chain, A.y);
          B.y = clampY(B.chain, B.y);
          moved++;
        }}
      }}
    }}
    if (moved === 0) break;
  }}

  for (const [rid, r] of Object.entries(RES)) {{
    const g = document.getElementById('res_' + slug(rid));
    if (!g) continue;
    g.setAttribute('transform', `translate(${{r.x}},${{r.y}})`);
  }}
}}

let view = {{x:0, y:0, k:1}};
function applyView() {{
  viewport.setAttribute('transform', `translate(${{view.x}},${{view.y}}) scale(${{view.k}})`);
}}

function svgPoint(evt) {{
  const pt = svg.createSVGPoint();
  pt.x = evt.clientX; pt.y = evt.clientY;
  const ctm = svg.getScreenCTM();
  if (!ctm) return {{x:0,y:0}};
  return pt.matrixTransform(ctm.inverse());
}}

function getBBoxVisible() {{
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const [rid, r] of Object.entries(RES)) {{
    if (hbOnly && !HB_RES.has(rid)) continue;
    const rad = nodeRadius(r) + 12;
    minx = Math.min(minx, r.x - rad);
    miny = Math.min(miny, r.y - rad);
    maxx = Math.max(maxx, r.x + rad);
    maxy = Math.max(maxy, r.y + rad);
  }}
  if (!isFinite(minx)) return {{minx:0,miny:0,maxx:W,maxy:H}};
  return {{minx, miny, maxx, maxy}};
}}

function fitViewToResidues() {{
  const bb = getBBoxVisible();
  const pad = 40;
  const bw = (bb.maxx - bb.minx) + 2*pad;
  const bh = (bb.maxy - bb.miny) + 2*pad;

  const k = Math.min(W / bw, H / bh);
  view.k = Math.max(0.2, Math.min(8.0, k));

  const cx = (bb.minx + bb.maxx) / 2;
  const cy = (bb.miny + bb.maxy) / 2;
  view.x = W/2 - cx * view.k;
  view.y = H/2 - cy * view.k;
  applyView();
}}

let draggingNode = null;
let draggingPan = null;
let rotating = null;

svg.addEventListener('pointerdown', (evt) => {{
  const p = svgPoint(evt);
  const g = evt.target.closest('g.res');

  const isMiddle = (evt.button === 1);
  const isShift = evt.shiftKey;

  if (panZoomEnabled && (isMiddle || isShift) && !g) {{
    draggingPan = {{x0: p.x, y0: p.y, vx0: view.x, vy0: view.y}};
    svg.setPointerCapture(evt.pointerId);
    return;
  }}

  if (!dragEnabled) return;
  if (!g) return;

  const rid = g.getAttribute('data-res');
  if (!rid || !RES[rid]) return;

  // if filtered out, do nothing
  if (hbOnly && !HB_RES.has(rid)) return;

  if (RES[rid].focus && evt.altKey) {{
    const sx = (p.x - view.x) / view.k;
    const sy = (p.y - view.y) / view.k;

    const r = RES[rid];
    const vx = sx - r.x;
    const vy = sy - r.y;
    const ang0 = Math.atan2(vy, vx);

    rotating = {{ rid, ang0, rot0: (r.rot || 0) }};
    g.setPointerCapture(evt.pointerId);
    return;
  }}

  const sx = (p.x - view.x) / view.k;
  const sy = (p.y - view.y) / view.k;
  draggingNode = {{ rid, dx: sx - RES[rid].x, dy: sy - RES[rid].y }};
  g.setPointerCapture(evt.pointerId);
}});

svg.addEventListener('pointermove', (evt) => {{
  const p = svgPoint(evt);

  if (draggingPan) {{
    const dx = (p.x - draggingPan.x0);
    const dy = (p.y - draggingPan.y0);
    view.x = draggingPan.vx0 + dx;
    view.y = draggingPan.vy0 + dy;
    applyView();
    return;
  }}

  if (rotating) {{
    const rid = rotating.rid;
    const r = RES[rid];
    if (!r) return;

    const sx = (p.x - view.x) / view.k;
    const sy = (p.y - view.y) / view.k;

    const vx = sx - r.x;
    const vy = sy - r.y;
    const ang = Math.atan2(vy, vx);

    const dAng = ang - rotating.ang0;
    r.rot = rotating.rot0 + dAng * 180.0 / Math.PI;

    setAtomBlockRotation(rid);
    updateAllEdges();
    return;
  }}

  if (!draggingNode) return;
  const r = RES[draggingNode.rid];
  if (!r) return;

  const sx = (p.x - view.x) / view.k;
  const sy = (p.y - view.y) / view.k;

  const nx = sx - draggingNode.dx;
  let ny = sy - draggingNode.dy;
  ny = clampY(r.chain, ny);

  r.x = nx; r.y = ny;

  const g = document.getElementById('res_' + slug(draggingNode.rid));
  if (g) g.setAttribute('transform', `translate(${{nx}},${{ny}})`);
  updateAllEdges();
}});

svg.addEventListener('pointerup', () => {{
  rotating = null;
  draggingNode = null;
  draggingPan = null;
  if (RELAX_ON_POINTERUP) {{
    resolveOverlaps();
    updateAllEdges();
  }}
}});
svg.addEventListener('pointercancel', () => {{
  rotating = null;
  draggingNode = null;
  draggingPan = null;
}});

svg.addEventListener('wheel', (evt) => {{
  if (!panZoomEnabled) return;
  evt.preventDefault();

  const p = svgPoint(evt);
  const zoomFactor = Math.exp(-evt.deltaY * 0.0016);

  const k0 = view.k;
  const k1 = Math.max(0.2, Math.min(8.0, k0 * zoomFactor));
  if (k1 === k0) return;

  const sx = (p.x - view.x) / k0;
  const sy = (p.y - view.y) / k0;

  view.k = k1;
  view.x = p.x - sx * k1;
  view.y = p.y - sy * k1;
  applyView();
}}, {{passive:false}});

function applyHbondVisibilityByFilter() {{
  for (const e of EDGES) {{
    if (e.type !== 'hbond') continue;
    const lid = 'h_' + slug(e.ra) + '_' + e.a + '__' + slug(e.rb) + '_' + e.b;
    const line = document.getElementById(lid);
    const text = document.getElementById(lid + '_t');
    const show = (!hbOnly) || (HB_RES.has(e.ra) && HB_RES.has(e.rb));
    if (line) line.style.display = show ? '' : 'none';
    if (text) text.style.display = show ? '' : 'none';
  }}
}}

resolveOverlaps();
for (const rid of Object.keys(RES)) setAtomBlockRotation(rid);
updateAllEdges();
fitViewToResidues();
applyView();
applyContactsVisibility();
applyResidueFilter();
applyHbondVisibilityByFilter();
{js_extra}
applyTheme();
applyAAMode();
// Initial focus label placement (avoid atom-label overlap on first load)
requestAnimationFrame(() => {{
  requestAnimationFrame(() => {{
    try {{ autoPlaceFocusTopLabelsAndPull(); }} catch (e) {{ console.warn('init focus label placement failed', e); }}
  }});
}});

</script>
</body>
</html>
"""

    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("drw", help="LigPlot+ .drw file")
    ap.add_argument("--out", default="2DPPIViewer.html", help="Output HTML")

    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--height", type=int, default=900)
    ap.add_argument("--margin", type=int, default=80)

    ap.add_argument("--atom_scale", type=float, default=3.0)
    ap.add_argument("--split_extra", type=float, default=35.0)

    ap.add_argument("--x_expand", type=float, default=1.50)
    ap.add_argument("--y_expand", type=float, default=1.00)
    ap.add_argument("--x_stretch", type=float, default=None)
    ap.add_argument("--y_stretch", type=float, default=None)

    ap.add_argument("--relax_iters", type=int, default=320)
    ap.add_argument("--relax_gap_px", type=float, default=14.0)
    ap.add_argument("--relax_step", type=float, default=0.38)
    ap.add_argument("--relax_on_pointerup", action="store_true")

    ap.add_argument("--simple_r", type=float, default=20.0)
    ap.add_argument("--focus_r", type=float, default=44.0)

    ap.add_argument("--boundary_extend_px", type=float, default=500.0)

    ap.add_argument("--atom_label_mode", choices=["element", "atomname"], default="atomname")
    ap.add_argument("--show_atom_labels_for_all", action="store_true")

    ap.add_argument("--chain_label_offset_px", type=float, default=26.0)

    args = ap.parse_args()

    if args.x_stretch is not None:
        args.x_expand = args.x_stretch
    if args.y_stretch is not None:
        args.y_expand = args.y_stretch

    build_html(
        drw_path=args.drw,
        out_html=args.out,
        width=args.width,
        height=args.height,
        margin=args.margin,
        atom_scale=args.atom_scale,
        split_extra=args.split_extra,
        x_expand=args.x_expand,
        y_expand=args.y_expand,
        relax_iters=args.relax_iters,
        relax_gap_px=args.relax_gap_px,
        relax_step=args.relax_step,
        simple_r=args.simple_r,
        focus_r=args.focus_r,
        relax_on_pointerup=args.relax_on_pointerup,
        boundary_extend_px=args.boundary_extend_px,
        atom_label_mode=args.atom_label_mode,
        show_atom_labels_for_all=args.show_atom_labels_for_all,
        chain_label_offset_px=args.chain_label_offset_px,
    )
    print(f"[OK] saved: {args.out}")


if __name__ == "__main__":
    main()
