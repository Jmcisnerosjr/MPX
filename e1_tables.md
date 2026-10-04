# E1 orbital calculations

Method: independent fixed Kepler ellipses from orbital_elements.json / network_model.json. M(t) = M0 + n t (days) mod 360; Kepler's equation by Newton to |dE| < 1e-15 rad; orbit-plane (x', y') rotated by omega, i, Omega into the Sun-centred J2000 ecliptic frame (AU). Time: hours since 2026-09-22 00:00 TDB, 86,400 s days. Light time solved by fixed-point iteration on t_a = t_e + 8.317 min/AU * |r_recv(t_a) - r_send(t_e)| until successive iterates differ by < 1e-9 h (3.6 us), 280x tighter than the brief's 1 ms.

## Epoch check (t = 0)

| Body | x | y | z | max abs error vs brief (AU) | Pass (<1e-5) |
|---|---|---|---|---|---|
| Mercury | -0.213426 | -0.410281 | -0.013955 | 1.57e-07 | yes |
| Venus | 0.679067 | -0.257950 | -0.042725 | 2.15e-07 | yes |
| Earth | 1.003581 | -0.023693 | -0.000003 | 3.84e-07 | yes |
| Mars | 0.247476 | 1.526869 | 0.025929 | 4.68e-07 | yes |
| Ceres | 0.375786 | 2.653941 | 0.014791 | 3.89e-07 | yes |
| Jupiter | -3.438179 | 4.038309 | 0.060149 | 3.98e-07 | yes |
| Saturn | 9.271221 | 1.717485 | -0.398962 | 3.39e-07 | yes |
| Uranus | 8.962475 | 17.253919 | -0.052129 | 4.51e-07 | yes |
| Neptune | 29.839098 | 1.351785 | -0.715470 | 2.75e-07 | yes |
| Relay A | 2.000000 | 2.000000 | 0.000000 | 4.44e-16 | yes |
| Relay B | -2.000000 | 2.000000 | 0.000000 | 4.44e-16 | yes |

Worst error: 4.68e-07 AU (limit 1e-5). Errors are at the 6-decimal rounding of the brief's table.

## Period return and radius bounds

| Body | Period (d) | abs(r(P) - r(0)) (AU) | r min sampled | a(1-e) | r max sampled | a(1+e) |
|---|---|---|---|---|---|---|
| Mercury | 87.9688 | 3.6e-16 | 0.307495 | 0.307495 | 0.466700 | 0.466700 |
| Venus | 224.7032 | 3.5e-16 | 0.718443 | 0.718443 | 0.728232 | 0.728232 |
| Earth | 364.9851 | 0.0e+00 | 0.982458 | 0.982458 | 1.016552 | 1.016552 |
| Mars | 686.9537 | 1.6e-15 | 1.381314 | 1.381314 | 1.665992 | 1.665992 |
| Ceres | 1680.2174 | 2.5e-15 | 2.545352 | 2.545352 | 2.986553 | 2.986553 |
| Jupiter | 4332.9562 | 7.5e-15 | 4.952740 | 4.952740 | 5.453457 | 5.453457 |
| Saturn | 10761.9640 | 2.6e-14 | 9.014072 | 9.014072 | 10.066883 | 10.066883 |
| Uranus | 30881.2692 | 8.8e-14 | 18.348367 | 18.348367 | 20.178834 | 20.178834 |
| Neptune | 60208.5952 | 0.0e+00 | 29.770212 | 29.770212 | 30.357890 | 30.357890 |
| Relay A | 1737.4644 | 2.7e-15 | 2.828427 | 2.828427 | 2.828427 | 2.828427 |
| Relay B | 1737.4644 | 2.7e-15 | 2.828427 | 2.828427 | 2.828427 | 2.828427 |

## Worked positions at later times (AU)

| Body | h 300 | +1 yr | +10 yr | +100 yr |
|---|---|---|---|---|
| Mercury | (0.05684, -0.45413, -0.04233) | (0.07625, -0.44929, -0.04371) | (0.17882, 0.25506, 0.00444) | (0.17390, -0.40655, -0.04917) |
| Venus | (0.72531, -0.01204, -0.04201) | (-0.65071, -0.31094, 0.03327) | (0.23911, 0.68011, -0.00445) | (-0.71792, 0.03327, 0.04188) |
| Earth | (0.98216, 0.18935, -0.00002) | (1.00360, -0.01915, -0.00000) | (1.00287, 0.02169, -0.00001) | (0.90472, 0.41692, -0.00004) |
| Mars | (0.08025, 1.56020, 0.03073) | (-0.48538, -1.41726, -0.01780) | (-1.65116, 0.00559, 0.04060) | (-1.17694, 1.16222, 0.05321) |
| Ceres | (0.24456, 2.65922, 0.03913) | (-2.48285, 0.37527, 0.46927) | (-2.22657, 1.14969, 0.44657) | (2.89975, -0.35361, -0.54538) |
| Jupiter | (-3.51062, 3.98088, 0.06201) | (-5.04394, 1.93586, 0.10480) | (1.43706, 4.85850, -0.05232) | (0.41756, -5.21621, 0.01231) |
| Saturn | (9.25439, 1.78584, -0.39949) | (8.57139, 3.65962, -0.40505) | (-7.40833, 5.44258, 0.19972) | (-9.00531, 2.50723, 0.31456) |
| Uranus | (8.91838, 17.27421, -0.05148) | (7.65190, 17.80133, -0.03316) | (-5.11057, 18.07449, 0.13284) | (-12.03444, 14.05995, 0.20753) |
| Neptune | (29.83697, 1.39118, -0.71624) | (29.75559, 2.50168, -0.73739) | (27.05294, 12.48151, -0.88206) | (-22.51059, -20.33742, 0.94031) |
| Relay A | (1.90758, 2.08833, 0.00000) | (-1.44315, 2.43255, 0.00000) | (0.40378, 2.79946, 0.00000) | (1.70515, 2.25665, 0.00000) |
| Relay B | (-2.08833, 1.90758, 0.00000) | (-2.43255, -1.44315, -0.00000) | (-2.79946, 0.40378, 0.00000) | (-2.25665, 1.70515, 0.00000) |

## Receiver motion changes arrival time

**Neptune -> Relay A, t_e = 1 s:** frozen-receiver arrival 3h51m41.618s, moving-receiver arrival 3h51m42.212s, difference **+0.594 s** (594x the 1 ms tolerance).

Fixed-point iterates (s after t_e): 13900.618008, 13901.212022, 13901.212047, 13901.212047, 13901.212047

**Earth -> Relay A, t_e = 1 s:** frozen-receiver arrival 0h18m46.639s, moving-receiver arrival 0h18m46.661s, difference **+0.021 s** (21x the 1 ms tolerance).

Fixed-point iterates (s after t_e): 1125.639333, 1125.660747, 1125.660748, 1125.660748, 1125.660748

**Mercury -> Relay B, t_e = 1 s:** frozen-receiver arrival 0h24m58.184s, moving-receiver arrival 0h24m58.171s, difference **-0.013 s** (13x the 1 ms tolerance).

Fixed-point iterates (s after t_e): 1497.184481, 1497.171480, 1497.171480, 1497.171480, 1497.171480

## Forward and return analysed separately

Round trip Earth gateway -> Relay A -> Earth at t_ready = 0: forward packet, relay 1 s processing + 1 s serialization, return launch from the relay's arrival-time position.

| Leg | t_e | t_a | flight (s) | d (AU) | Sun clearance (AU) |
|---|---|---|---|---|---|
| Earth -> Relay A | 0h00m01.000s | 0h18m46.661s | 1125.661 | 2.255743 | 1.0039 |
| Relay A -> Earth | 0h18m48.661s | 0h37m34.121s | 1125.460 | 2.255341 | 1.0039 |

Forward and return flights differ by 0.201 s because both endpoints move during the round trip; a symmetric-path assumption would be wrong by this much.
