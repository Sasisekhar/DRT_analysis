#!/usr/bin/env python3
"""Draw DRT task graphs and demand bound functions from `drt --json` output.

Usage:
    ./drt --json 60 | python3 draw_drt.py            # open a window per figure
    ./drt --json 60 > tasks.json
    python3 draw_drt.py tasks.json -o figures        # write figures/*.png

Requires matplotlib. networkx (+ scipy) is optional and gives nicer layouts;
without it vertices are placed on a circle.
"""
import argparse
import json
import math
import os
import re
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.path import Path

R = 0.11  # vertex radius in layout units


# ------------------------------------------------------------------ geometry

def unit(x, y):
    n = math.hypot(x, y)
    return (x / n, y / n) if n > 0 else (0.0, 1.0)


def rotate(v, deg):
    a = math.radians(deg)
    return (v[0] * math.cos(a) - v[1] * math.sin(a), v[0] * math.sin(a) + v[1] * math.cos(a))


def add(p, v, s=1.0):
    return (p[0] + s * v[0], p[1] + s * v[1])


def layout(n, edges):
    if n == 1:
        return {0: (0.0, 0.0)}
    try:
        import networkx as nx
        g = nx.Graph()
        g.add_nodes_from(range(n))
        g.add_edges_from((e["from"], e["to"]) for e in edges if e["from"] != e["to"])
        pos = nx.kamada_kawai_layout(g)
        return {i: (float(pos[i][0]), float(pos[i][1])) for i in range(n)}
    except Exception:  # networkx/scipy missing: fall back to a circle
        return {i: (math.cos(2 * math.pi * i / n + math.pi / 2),
                    math.sin(2 * math.pi * i / n + math.pi / 2)) for i in range(n)}


# ------------------------------------------------------------------ task graph

def edge_label(ax, xy, text):
    ax.text(*xy, text, ha="center", va="center", fontsize=9, zorder=4,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))


def draw_edge(ax, pu, pv, rad, p):
    # arc3 places its control point at mid + rad * (dy, -dx); use the tangents
    # towards it to start/end the arrow on the vertex circles.
    dx, dy = pv[0] - pu[0], pv[1] - pu[1]
    a = add(pu, unit(dx / 2 + rad * dy, dy / 2 - rad * dx), R)
    b = add(pv, unit(dx / 2 - rad * dy, dy / 2 + rad * dx), -R)
    ax.add_patch(FancyArrowPatch(a, b, connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
                                 mutation_scale=14, lw=1.2, color="0.25", shrinkA=0, shrinkB=0,
                                 zorder=2))
    dx, dy = b[0] - a[0], b[1] - a[1]
    mid = ((a[0] + b[0]) / 2 + 0.5 * rad * dy, (a[1] + b[1]) / 2 - 0.5 * rad * dx)
    edge_label(ax, mid, str(p))


def draw_self_loop(ax, pv, out, p):
    s, t = rotate(out, 35), rotate(out, -35)
    path = Path([add(pv, s, R), add(pv, rotate(out, 40), 3.4 * R),
                 add(pv, rotate(out, -40), 3.4 * R), add(pv, t, R)],
                [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=14, lw=1.2,
                                 color="0.25", zorder=2))
    edge_label(ax, add(pv, out, 2.9 * R), str(p))


def draw_graph(ax, task):
    V, E = task["vertices"], task["edges"]
    pos = layout(len(V), E)
    cx = sum(x for x, _ in pos.values()) / len(pos)
    cy = sum(y for _, y in pos.values()) / len(pos)
    outward = {i: unit(x - cx, y - cy) for i, (x, y) in pos.items()}

    directed = {(e["from"], e["to"]) for e in E}
    seen = {}
    loops = set()
    for e in E:
        u, v = e["from"], e["to"]
        k = seen.get((u, v), 0)  # parallel edges in the same direction fan out
        seen[(u, v)] = k + 1
        if u == v:
            loops.add(u)
            draw_self_loop(ax, pos[u], rotate(outward[u], 50 * k), e["p"])
            continue
        rad = (0.2 if (v, u) in directed else 0.0) + 0.2 * k
        draw_edge(ax, pos[u], pos[v], rad, e["p"])

    for i, v in enumerate(V):
        ax.add_patch(Circle(pos[i], R, fc="#e8eef8", ec="#3b5b92", lw=1.5, zorder=3))
        ax.text(*pos[i], v["id"], ha="center", va="center", fontsize=10, zorder=5)
        side = rotate(outward[i], 90) if i in loops else outward[i]
        ax.text(*add(pos[i], side, 2.1 * R), f"⟨{v['e']}, {v['d']}⟩", ha="center",
                va="center", fontsize=9, color="#3b5b92", zorder=5)

    xs = [x for x, _ in pos.values()]
    ys = [y for _, y in pos.values()]
    m = 5 * R
    ax.set_xlim(min(xs) - m, max(xs) + m)
    ax.set_ylim(min(ys) - m, max(ys) + m)
    ax.set_aspect("equal")
    ax.axis("off")
    u = task["utilization"]
    kind = "frame-separated" if task["frame_separated"] else "arbitrary deadlines"
    ax.set_title(f"{task['name']}   U(T) = {u['e']}/{u['p']}   ({kind})")


# ------------------------------------------------------------------ dbf plots

def staircase(steps, horizon):
    xs, ys = [0], [0]
    for d, e in steps:
        xs.append(d)
        ys.append(e)
    xs.append(horizon)
    ys.append(ys[-1])
    return xs, ys


def draw_dbf(ax, task):
    steps, horizon = task["dbf"]["steps"], task["dbf"]["horizon"]
    xs, ys = staircase(steps, horizon)
    ax.step(xs, ys, where="post", color="#3b5b92", lw=1.5)
    if steps:
        ax.plot([d for d, _ in steps], [e for _, e in steps], "o", color="#3b5b92", ms=4)
    ax.set_xlim(0, horizon)
    ax.set_ylim(0, max(ys) * 1.15 + 1)
    ax.set_xlabel("t")
    ax.set_ylabel("dbf_T(t)")
    ax.set_title("demand bound function")
    ax.grid(alpha=0.3)


def draw_task_set(ax, data):
    tasks, fz = data["tasks"], data["feasibility"]
    horizon = max(t["dbf"]["horizon"] for t in tasks) if tasks else 1
    events = {}
    for t in tasks:
        prev = 0
        for d, e in t["dbf"]["steps"]:
            events[d] = events.get(d, 0) + e - prev
            prev = e
    xs, ys, total = [0], [0], 0
    for d in sorted(events):
        total += events[d]
        xs.append(d)
        ys.append(total)
    xs.append(horizon)
    ys.append(total)

    ax.step(xs, ys, where="post", color="#3b5b92", lw=1.5, label="dbf(t)")
    ax.plot([0, horizon], [0, horizon], "--", color="0.5", lw=1, label="t")
    if fz["bound"] > 0:
        ax.axvline(fz["bound"], color="#c77c02", ls=":", label=f"bound = {fz['bound']}")
    if fz["witness"] is not None:
        w = fz["witness"]
        ax.plot([w], [ys[xs.index(w)]], "o", color="#c0392b", ms=8, label=f"witness t_f = {w}")
    ax.set_xlim(0, horizon)
    ax.set_ylim(0, max(max(ys), horizon) * 1.05)
    ax.set_xlabel("t")
    ax.set_ylabel("demand")
    ax.set_title(f"task set: {fz['status']}   U(τ) = {fz['utilization']:.4g}\n{fz['reason']}",
                 fontsize=10)
    ax.legend(loc="upper left")
    ax.grid(alpha=0.3)


# ------------------------------------------------------------------ main

def slug(name, i):
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", name).strip("_")
    return s or f"task{i}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("json", nargs="?", help="output of `drt --json` (default: stdin)")
    ap.add_argument("-o", "--out", help="directory to write PNGs to instead of showing windows")
    args = ap.parse_args()

    data = json.load(open(args.json) if args.json else sys.stdin)
    figures = []
    for i, task in enumerate(data["tasks"]):
        fig, (g, d) = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.2, 1]})
        draw_graph(g, task)
        draw_dbf(d, task)
        fig.tight_layout()
        figures.append((slug(task["name"], i), fig))
    fig, ax = plt.subplots(figsize=(8, 5.5))
    draw_task_set(ax, data)
    fig.tight_layout()
    figures.append(("task_set", fig))

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        for name, fig in figures:
            path = os.path.join(args.out, f"{name}.png")
            fig.savefig(path, dpi=150)
            print(path)
    else:
        plt.show()


if __name__ == "__main__":
    main()
