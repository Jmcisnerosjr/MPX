# MPX: The MultiPlanetary Exchange

A stock exchange for traders spread across the solar system, where every message moves at the speed of light. From Neptune, that's over four hours one way, and almost half the messages get lost.

**Gator Quant Hacks 2026, Quantitative Puzzles track.** Jose Cisneros and Benjamin Robertson.

**Live app:** https://claude.ai/artifact/2kxopKxCWDEvPpsZcSVgLY

## What's here

| Folder | Contents |
| --- | --- |
| `paper/` | The submission: `MPX_Cisneros_Robertson.pdf` (12-page design paper + 12-page evidence appendix) and its LaTeX source |
| `app/` | `index.html`, the interactive app. Open it in any browser; it runs offline |
| `toolkit/` | Python model of the brief, the 200-year scan, the protocol engine, stress tests and figures |
| `toolkit/data/` | The brief's frozen v1.3 data (orbital elements, network model, JPL Horizons rows) |
| `toolkit/out/` | Precomputed results, including `all_messages.csv`, every message of every scenario run |

## The design in brief

- **Two relays, never both blocked.** Relays A and B share one orbit 90 degrees apart, so the Sun can never block both from any planet. Checked every hour for 200 years.
- **Fair batches.** Orders clear at Earth in batches whose deadlines are set by the straight-line speed of light, so nobody trades on news that hadn't reached them.
- **Copies, not retries.** Deadline messages go out as several copies, so at least one arrives on time 99.5% of the time, within the 600-message daily limit.
- **Exactly-once money.** Only the destination can declare a transfer credited or never credited. Lost packets cost time, never money.
- **Margin for the slowest round trip,** with a fat-tail multiplier and a prefunded default waterfall.

## Run the app

Open `app/index.html` in a browser. Tabs: Home (portfolio), Market, Make a trade, Network, Risk.

## Reproduce the paper's numbers

Requires Python 3 with `numpy`, `scipy` and `matplotlib`.

```bash
cd toolkit
python3 verify.py      # 8 independent checks, all must PASS
python3 e1_orbits.py   # E1 orbits and light-time
python3 e4_scan.py     # E4 200-year closure scan (about 4 minutes)
python3 e4_analyze.py
python3 hub_loss.py
python3 s3_e2_e5.py    # S3 access tables, E2 probabilities, E5 shifted epochs
python3 design.py      # copies, deadlines, send-early times
python3 tiers.py       # two-speed schedule and daily quota over 200 years
python3 scenarios.py   # S1 traces
python3 runs_all.py    # S2 incident runs, S3 timelines, E5 reruns
python3 e5_runs.py     # E5 Monte Carlo price runs
python3 dvp_basis.py   # swaps, basis bands, expiries, receipts
python3 phase.py       # failure boundary grid
python3 figs.py        # paper figures
```

The model is also implemented independently in JavaScript (the app), and both agree to the last digit at hour 0: Neptune to Earth in 4.1747 h, 54.76% first-try success, 3.33% abandonment.

## AI use

We used Claude (Anthropic) as a coding and writing assistant throughout. Every number in the paper comes from the scripts in this repository.

## License

MIT
