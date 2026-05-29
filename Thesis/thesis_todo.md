# Thesis TODO

A working checklist of everything still outstanding, grouped by chapter.
Pulled from the inline `\todo{}` markers in the rewritten methodology and
results, plus things that are clearly missing or inconsistent across the
draft.

---

## High-priority (numbers / consistency / blocking)

- [ ] **Reconcile reeds TPR/FPR numbers.** Prose says March TPR = 0.93,
      April TPR = 0.57 at FPR = 0.67. Table 6 says March CR1 TPR = 0.88,
      April CR1 TPR = 0.65, FPR ~0.6. Pick one set and use it everywhere
      (text, table, figure captions).
- [ ] **Reconcile control-scan clutter floor numbers.** Prose says
      $-27.8$\,dB, table says Own clutter floor $= -36.9$\,dB, and the
      fixed detection threshold is $-40.1$\,dB. The prose is
      self-consistent (66 - (-27.8) = 93.8 dB matches CR1 SNR), but three
      different floor values are floating around. Settle on definitions
      (clutter floor vs threshold vs table value) and make them match.
- [ ] **Pull chirp bandwidth $B$ and sweep time from the SiRad
      controller config.** Chip-level specs (TX power, NF, VCO tuning
      range, antenna gain) are now filled in from the TRA\_120\_002
      datasheet. The chirp and sweep values are host-side and still
      need to go in so that $\Delta R = c/2B$ and the unambiguous range
      can be stated.
- [ ] **Read the $-3$\,dB beamwidth off Figs.\ 17 and 18 of the
      TRA\_120\_002 datasheet** and state it in the Radar System
      subsection. Then compute the module-level angular resolution
      using whatever lens/horn is in use.
- [ ] **Confirm whether a dielectric lens or horn was fitted** to the
      radar front end during these experiments. If yes, add its
      gain, focal length, and the resulting effective beamwidth; the
      $\sim$10\,m bare-chip range limit is otherwise right at the
      scans' range extent.
- [ ] **Convert control-scan attenuation slope into a directly comparable
      dB/m figure** so the Cherry Ballart 6.27 dB/m can be stated as a
      proper relative number ("X times higher than free-space").
- [ ] **Ground-truth check on the vegetation loss numbers**
      (March 58.3 dB, April 19.2 dB). Sanity-check the free-space
      propagation model used to derive expected amplitudes at 10.92 m
      and 13.88 m.

---

## Methodology

### Hardware
- [ ] Radar System: chirp bandwidth, sweep time, TX power, antenna type +
      model + beamwidth, ADC sample rate, manufacturer, scan rate.
- [ ] Derive and state $\Delta R$, max unambiguous range, angular
      resolution from the above.
- [ ] LiDAR: confirm the firmware/version of the Mid-70 used, since some
      of the spec numbers depend on it.
- [ ] GPS: state base-station ID/location, baseline length, and average
      satellite count during fixes.
- [ ] Sensor Mounting: add a photo or drawing of the mount with measured
      sensor offsets between radar and LiDAR. Even a 5 cm baseline matters
      at 5–10 m ranges.
- [ ] Corner Reflectors: tabulate physical dimensions, theoretical RCS
      at 120 GHz, theoretical LiDAR cross-section of the mounted mirror.

### Plane Tree
- [ ] Tree geometry: add an annotated figure showing the 15.5 × 18 × 20 m
      envelope.
- [ ] Scan procedure: tabulate scan positions, headings, dates, weather,
      ground conditions, tripod heights.
- [ ] Tree-specific preprocessing: write up the actual MATLAB pipeline
      values used (height filter bounds, $k$, multiplier, bounding-box
      dimensions, amplitude floor).
- [ ] Vertical alignment: revisit once tripod heights for the two
      sessions are recovered (or document that they're permanently lost).

### Reeds
- [ ] Table caption for `tab:scan_params` is referenced from methodology
      but the table itself lives in results — either move the table or
      duplicate it / add a forward reference.
- [ ] Behind-foliage region: add a small diagram showing which slice of
      the scan was used.
- [ ] Document why CFAR was rejected with a one-line numerical
      justification (sample density per range bin vs CFAR window
      requirements).

### Shrubs
- [ ] Cherry Ballart: write up the multi-azimuth cone analysis in full
      once it is run. Include sensitivity to cone half-angle.
- [ ] Define "cone half-angle" relative to the radar's actual beamwidth
      (i.e. is 15° narrower or wider than the beam?).
- [ ] LiDAR preprocessing: write the actual pipeline used (cropping,
      downsampling, intensity filtering, statistical outlier removal,
      ground removal). Currently a placeholder.

### Detection metrics
- [ ] Add the explicit equation for two-way attenuation estimate.
- [ ] State the free-space model derived from the control scan
      (slope + intercept in dB vs metres).
- [ ] Add a short "uncertainty budget" subsection — repeatability, beam
      divergence, angular sampling, GPS misalignment, weather. This is
      promised in the Aims and Objectives.

---

## Results

### Control scan
- [ ] Add a one-line statement of the free-space propagation model fit
      (slope + intercept + R²) so it can be referenced directly when the
      vegetation loss is computed later.
- [ ] Reconcile Reflector 1 and Reflector 2 LiDAR centroids
      ((8.71, -1.95, 0.33) and (9.66, 2.77, 0.42)) with the GPS-measured
      reflector positions, and show the two are consistent.

### Reeds
- [ ] Fix the TPR/FPR mismatch (see High-priority above).
- [ ] Add an explicit comment on the gravel multipath: ideally, run a
      quick sanity-check sim or a back-of-envelope calculation to confirm
      the +64 dB jump is plausibly multipath-driven rather than something
      else.
- [ ] Add a small figure caption for `fig:lidar_reflectivity` (currently
      "Enter Caption").
- [ ] Replace the duplicate `\label{fig:placeholder}` labels — there are
      at least three figures using the same label, which will silently
      break cross-references.
- [ ] Add a short paragraph explicitly stating that LiDAR fully fails to
      detect the reflector in the reeds, and tie this back to the
      research question.

### Cherry Ballart
- [ ] Write up the full results — currently mostly figure dumps and
      bullet points.
- [ ] Multi-azimuth cone sweep: show how attenuation varies with
      direction.
- [ ] Compare with the LiDAR result for the same shrub.
- [ ] Add a short paragraph on the fence multipath visible in
      `fig:visiblecr`.

### Oldschool Building bushes
- [ ] Currently a figure dump. Write up scan parameters, clutter floor,
      CR detection result, comparison with the Ballart and reeds.
- [ ] Decide whether this scan adds anything beyond Ballart — if not,
      move it to an appendix.

### 11/04/2026 shrub scans
- [ ] Currently only scan parameters. Add scene description,
      photographs, what shrub was scanned, results.
- [ ] State the actual beamwidth and confirm overlap = beamwidth / 2.

### Plane Tree (radar vs LiDAR)
- [ ] Add a quantitative comparison: per-bin radar/LiDAR return ratio as
      a function of range into the canopy, as the supervisor suggested
      in the conclusion notes.
- [ ] Compute and state radar attenuation through the canopy in dB/m so
      it can be compared with the Ballart and reeds numbers.

---

## Other chapters (not part of this rewrite, but flagged)

### Front matter
- [ ] Contribution of Work — empty.
- [ ] Abstract — empty.
- [ ] Acknowledgments — empty.

### Introduction
- [ ] "FIND THE SIGNAL IN THE NOISE" placeholder under Introduction —
      either turn into a real epigraph or remove.
- [ ] `\todo{put this in an image}` at end of Thesis Outline — either
      build the chapter-flow diagram or remove the marker.

### Literature Review
- [ ] All subsections except "LiDAR in Vegetation Sensing → Limitations"
      are still empty headings: Foliage Penetration Radar (history),
      Microwave Propagation, Transition to MMW, Scattering Regimes,
      Clutter, Comparative Studies, GPS and Spatial Registration,
      Identified Research Gap.
- [ ] `\textcolor{red}{TODO: Add somethign about LiDAR reflector}` in
      LiDAR Limitations.
- [ ] Decide whether the "complementary sensing behaviour has motivated
      multimodal fusion" sentence should sit in the literature review
      *and* in the methodology, or only one of them.

### Discussion
- [ ] Stub — write up.

### Conclusion
- [ ] Currently contains a red-text supervisor brief, not actual text.
      Move the brief into a working notes file and write a real
      conclusion.

### Bibliography
- [ ] Confirm every `\cite{}` resolves. There are several in the
      introduction (e.g. `9760104`, `pauli2017mmwave`,
      `iturp67613`, `livox_manual`, `dodRadarHandbook`) — verify each
      is in `bibliography.bib` and uses the correct key.
- [ ] Decide on a citation style and stick to it (some keys are
      lowercase, some are mixed — check IEEEtr is happy).

---

## Style / housekeeping

- [ ] Run a spell-check pass — there are recurring typos: "perofrm",
      "perofrmed", "perfomred", "anlsyis", "anaysis", "atteniated",
      "signficantly", "Refernce", "approximatly", "minim", "sysrem",
      "tst cloud", "olderschool", "Lidar"/"liDAR" inconsistency.
- [ ] Pick one capitalisation: "LiDAR" everywhere (not "Lidar" or
      "lidar" or "liDAR").
- [ ] Pick one of "corner reflector" / "CR" early and stick to it.
- [ ] Strip stray `\textcolor{red}{...}` working notes once the relevant
      content has been resolved.
- [ ] Resolve duplicate `\label{fig:placeholder}` figures — these will
      cause silent cross-reference bugs.
- [ ] Decide on date format: "11~March~2026" vs "11/03/2026" vs
      "8~April~2026" — currently mixed.
- [ ] Captions I had to invent during the rewrite (originals were "Enter
      Caption"): `fig:lidar_control_amplitude`, `fig:lidar_control_height`
      ("Coloured by amplitude" / "Coloured by height"), `fig:height` /
      `fig:reflectivity` ("Coloured by height" / "Coloured by
      reflectivity"), `fig:cr_higher` / `fig:seethrough`. Replace with
      what the figures actually show.
- [ ] Section numbering: `secnumdepth` is set to 0, so subsections are
      not numbered. Confirm this is intentional. If it is, the inline
      "Chapter 1 introduces…" prose in the Thesis Outline reads
      strangely.
