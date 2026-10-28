# GlueForge

## Start here

**Double-click `GlueForge.bat`** (Windows) or **`GlueForge.command`** (macOS),
or run `py glueforge.py`. A page opens in your browser and everything happens
there: what is installed, buttons to install what is not, and the screens for
uploading a pair and working through it.

Only Python has to be there first. `numpy` and `scipy` install from the
landing page, and between them they cover every stage, docking included —
the built-in FFT engine needs nothing else. Two optional pieces also install
from that page: **RDKit**, which gates the chemistry stage and nothing else,
and **LightDock**, a second docking engine. Neither is needed for a verdict,
so a failed install costs you one stage rather than the tool.

The two engines have not been compared on this benchmark — the built-in one
is measured (`cache/recovery`), LightDock is not — so the choice between
them is offered without a recommendation. LightDock has no Windows wheel at
any version; `setup-wsl.sh` is the route to it, and the landing page points
there rather than suggesting an install that cannot work.

**The fastest way to see it work** is the examples on the landing page: five
solved glue complexes with the anchor, the ligand and the partner already
assigned. Loading one skips straight to a verdict on a crystal structure,
which is what every threshold in the tool was calibrated on.

Docking and the verdict are two steps on purpose. Docking proposes a few
hundred placements of protein B; the gate judges one of them. The dock screen
has an **Assess this pose** button that runs the gate on whichever pose is
selected — that is the whole of the second step.

A verdict reached that way says so on screen, and it should be read
accordingly. On the six bound-bound cases in `cache/recovery`, the pose the
built-in engine ranks first was 2.6 Å from the crystal placement once (6H0F)
and 20–35 Å out the other five times, against a cavity that tolerates about
2 Å. Docking here is not yet good enough to hand the gate something worth
judging; the examples exist so the gate can be seen working on placements
that are right.

A finished run exports as a single HTML file that opens on a machine with no
Python and nothing installed. That is the copy you send to someone.

Everything below is for working with the pipeline directly.

---

A pipeline and dashboard for the anchor-first glue-discovery workflow. It runs
the workflow on solved glue systems and reports two things: whether frustration
actually marks where the glue sits, and whether the composite cavity at that
site is a usable design target.

The measurements are the deliverable. The dashboard is how you inspect them.

---

## What it found

Fourteen crystal ternary complexes, glue deleted, frustration computed across
the anchor–partner interface. Negative means contacts next to where the glue
sat are more frustrated than the rest of that interface — the direction the
hypothesis predicts.

Run on the **published AWSEM-MD (2015) tables**, contact and burial terms.

| System | structures | mean Δ | direction |
|---|---|---|---|
| CDK12–cyclin K × DDB1 | 8 | −0.197 | 8/8 |
| DCAF15 × RBM39 (Faust form) | 3 | −0.193 | 3/3 |
| DCAF15 × RBM39 (Bussiere form) | 1 | +0.103 | 0/1 |
| CRBN × C2H2 zinc finger | 2 | +0.038 | 1/2 |

12 of 14 in the predicted direction, sign test p = 0.0065.

The threshold-free version is the stronger evidence and it is the part that
held up best. The rank correlation between a contact's frustration and its
distance to the ligand is positive in **13 of 14** structures (mean Spearman
rho +0.164, sign test p = 0.0009), with 8 of 14 individually significant. So
frustration rises continuously toward the ligand rather than stepping at an
arbitrary line.

### What changed when the published tables replaced the stand-in

The project ran for a while on a transparent stand-in potential assembled from
hydropathy and charge, because the AWSEM file was not to hand. Swapping in the
real tables - and adding the burial term, which the stand-in did not have -
moved the numbers:

| | stand-in | AWSEM 2015 |
|---|---|---|
| split comparison, direction | 12/14, p = 0.0065 | 12/14, p = 0.0065 |
| split comparison, mean effect | −0.332 | **−0.141** |
| burial-stratified | 12/14, p = 0.0065 | **10/14, p = 0.090** |
| gradient, direction | 13/14, p = 0.0009 | 13/14, p = 0.0009 |
| gradient, mean rho | +0.227 | +0.164 |
| individually significant gradients | 11/14 | 8/14 |

Read that honestly: **the direction is robust to the potential, the magnitude
is roughly halved, and the burial control no longer clears p = 0.05.** The two
potentials agree at r = 0.868 on which systems show the effect, so they are
measuring the same thing with different sensitivity. Anyone quoting this
should quote the AWSEM column.

Use `--gamma builtin` on any script to reproduce the stand-in numbers and see
the dependence for yourself.

**Read it with four caveats.** The run above used the stand-in potential, not
the published AWSEM tables — see *Potentials* below. Fourteen structures are
not fourteen independent observations: the eight CDK12 entries are one
interface with eight ligands, so across biological systems the count is closer
to two of three. The two DCAF15 crystal forms disagree in sign on the same
biological pair, which bounds reproducibility more honestly than any p-value
here. And the CRBN interfaces are small and glue-dominated — 46 of 73 contacts
are proximal in 6H0F — so that null may be a power problem rather than a real
absence.

### The cavity is the better design target

Frustration has no ligand term, so a molecule can only move it by changing
where the two proteins sit. The cavity at the glued interface does not have
that problem: a molecule placed there is scored directly.

The detector finds it. 11 of 14 systems recover 100% of the real ligand's
heavy atoms inside the detected cavity. The three failures are all the DCAF15
Faust form, whose enclosure quality is 0.05-0.25 - a surface groove rather than
a pocket, which is what the paper describes. Worth knowing, because
pocket-conditioned generative models are trained on enclosed pockets and a
groove is out of distribution for them.

The error budget, from displacing the partner rigidly:

Averaged over all 14 systems:

| partner error | shape (Jaccard) | r(electrostatic) | anchor side kept | partner side kept |
|---|---|---|---|---|
| 0.5 A | 0.86 | 0.98 | 0.92 | 0.89 |
| 1 A | 0.83 | 0.95 | 0.94 | 0.81 |
| 2 A | 0.72 | 0.91 | 0.89 | 0.55 |
| 3 A | 0.57 | 0.82 | 0.82 | 0.37 |
| 5 A | 0.51 | 0.76 | 0.80 | 0.36 |

Fields are far more robust than shape, because they are dominated by the
anchor side and by long-range terms. Design consequence: condition on the
field and weight shape loosely, which is the opposite of the default.

And the cavity is asymmetric. **61-86% of it is anchor-lined in every system**,
and that part comes from a co-crystal rather than from docking. The last two
columns above are the point: by 3 A the anchor-lined region still holds 0.82
while the partner-lined region has collapsed to 0.37. Most of the design target
is exact, and only the smaller rest carries your docking error.

### One convention for hydrogens, because the depositors do not agree

Nine of the fourteen structures model hydrogens and five do not. That is a
choice made by whoever deposited each one, and it was silently changing the
answer: a ligand's own hydrogens shield its own heavy atoms enough to move its
apparent exposed area by **1.7x to 4.3x**. The presented-face column was
therefore not comparable across the benchmark, and the same flaw sat in cavity
detection, where a modelled hydrogen projects about 1.1 A past its parent atom
and shrinks every cavity it touches.

It is the ligand's own hydrogens that do this. Dropping the protein's moves
the answer by about 2%. Radius conventions do not rescue it either — ProtOr,
NACCESS and Bondi agree within 3% — and a protein-sized probe makes it worse
rather than better, because the with-hydrogen measurement collapses toward
zero faster than the without.

So hydrogens are discarded everywhere, in `atom_sasa` and `detect_cavity`
alike. The cost is a little physical accuracy on the structures that have
them; the gain is that two structures can be compared at all.

What this corrected, for 6TD3:

| | before | after |
|---|---|---|
| presented area | 58.8 A^2 | 164.1 A^2 |
| buried fraction | 0.764 | 0.749 |

The ranking moved more than the magnitudes did. 6TD3 and 8BUJ used to read as
the most buried ligands in the set, at 59 and 21 A^2; they were simply the
ones with hydrogens. Under one convention 6Q0R/6Q0V/6Q0W separate out as
genuinely exposed (36-42% buried) and the rest sit at 70-83%.

Two results were checked rather than assumed, and both survive. **Frustration**
never saw a hydrogen — it runs on one representative atom per residue — so the
12/14 split and the 13/14 gradient are untouched. **Cavity recovery** is still
11/14 at 100%, and the three misses are exactly the three shallow-groove
ligands, which is what it should be: an unenclosed ligand has no enclosed
cavity around it.

Prefer **buried fraction** when comparing systems. It is stable across both
conventions to about +/-0.05, because the shielding largely cancels between
the bound and free states. Absolute area is not.

### Does the presented-face restraint actually carry signal?

This is the question the whole method rests on, and it turns out not to need
a docking program. Take each solved complex, throw the partner away, then put
it back in 400 uniform random orientations, each pushed in along a random
direction until the two just touch. Every decoy is a plausible complex --
contacting, not interpenetrating, not floating in space -- which is the only
kind worth comparing against. Score each on how much of the presented face it
buries and ask where the crystal pose lands.

`python3 sidecar/run_decoys.py`

The native pose ranks first in 13 of 14 systems and second in the other.
Median decoy burial is 0.0% in every system: most random placements, however
plausible as complexes, do not touch the presented face at all. In 12 of 14,
no decoy clears 50% burial at all.

Run at two decoy-contact distances, because that parameter had already turned
out to matter once. At 2.6 A -- how close the real partners come in these
crystals -- and at 3.0 A, van der Waals contact, the result is the same: 13 of
14 rank first, 14 of 14 in the top five. The conclusion does not depend on
the choice.

Three details of how the decoys are made turned out to matter more than the
result did, and two of them were wrong at first.

**Decoys have to be pushed in to the same distance the real partner sits at.**
The first version stopped at 1.5 A between heavy-atom centres, which is
interpenetration rather than contact, and let decoys bury flat faces they
could not really reach. Under that version the three groove-bound ligands
looked like failures -- 6Q0R ranked 9th, 6Q0W 19th. At proper contact they
rank 1st and 1st. The real partners come to 2.41-2.75 A in these crystals, so
`--gap` defaults to 2.6: holding decoys further out than the native sits
would hand the native burial it has not earned.

**Orientations have to be uniform over SO(3).** Three independent Euler
angles concentrate near the poles and would over-sample one band of the
anchor's surface. Sampled from normalised quaternions instead, and tested.

**A pose in the ligand's own pocket is not a result.** LightDock is handed
protein only -- its scoring function has no parameters for an arbitrary small
molecule -- so the pocket is empty as far as the docking run knows. A pose
that settles into it buries the presented face perfectly by occupying it, and
ranks first in an unfiltered run. Every pose is re-scored against the ligand.

### Picking restraints by distance to the ligand does not work

The obvious rule is that the active restraints should be the rim residues
nearest the exposed ligand. It is wrong, and it fails in the worst possible
way: it produces an active set the *correct answer* cannot satisfy, and
LightDock's own filter then discards the right pose.

Measured across all 324 rim residues in the benchmark, only 30-57% of them
are touched by the real partner, and which ones are is not predicted by
distance to the ligand at all:

| predictor | pooled r | p |
|---|---|---|
| distance to nearest exposed ligand atom | −0.07 | 0.22 |
| residue solvent accessibility | +0.33 | <1e-4 |
| distance from the face centroid | −0.40 | <1e-4 |
| **projection along the face's outward normal** | **+0.57** | **<1e-4** |

The reason is geometric. A partner approaches along the face normal, so two
residues equally close to the ligand are not equally reachable -- one is on
the approach side and the other is behind it. Ranking by normal projection
and taking the top six raises the crystal pose's own restraint satisfaction
from a median of 0.52 to 1.00. Under the old rule 6TD3's crystal pose
satisfied 0.38 of its restraints, just below LightDock's default 0.4 cut: the
right answer would have been filtered out.

### A face that wraps the ligand has no approach direction

Ranking by the normal only means something if there is a normal. Measuring
how much the exposed face points one way -- the area-weighted mean of the
unit vectors from the ligand centroid to each exposed atom -- separates the
benchmark cleanly:

| | directionality | |
|---|---|---|
| groove-bound (6Q0R, 6Q0V, 6Q0W) | 0.215 – 0.274 | 36-42% buried |
| pocket-bound (the other eleven) | 0.343 – 0.877 | 68-83% buried |

Below 0.30 the tool sets **no active restraints at all** and submits the whole
rim as permitted. The region still prunes the search; nothing unsatisfiable is
demanded of the right answer. Forcing an active set on a groove face is how
6Q0W ended up with a crystal pose satisfying none of its own restraints.

The threshold sits in an observed gap, but the gap is n=14. It is a soft
warning that changes how the restraint is written, not a verdict on the pair.
Directionality correlates with buried fraction at rho = +0.84, so the upload
screen can say which case you are in before anything is docked.

### The electrostatics were reporting the protein, not the site

The cavity's electrostatic term is a Coulomb sum over formal charges with a
distance-dependent dielectric. It had no distance cutoff, so it summed every
charge in both proteins -- which means it reported the proteins' overall
acidity as though it were a property of the site.

How much of it was not local:

| | no cutoff | ≤15 Å | ≤12 Å | ≤10 Å |
|---|---|---|---|---|
| 6TD3 | −3.45 | −0.94 | −0.47 | −0.18 |
| 8BU7 | −2.99 | −0.81 | −0.16 | +0.29 |
| 6H0F | **+0.75** | −0.70 | **−0.88** | −0.80 |

At 6TD3, 86% of the field came from charges more than 12 Å away. At 6H0F the
unbounded sum had the opposite sign to the local environment. Real solvent and
counterions screen long-range fields almost completely, and a
distance-dependent dielectric only partly stands in for that, so the sum is
now cut off at 12 Å.

This was found by following a suggestion that looked wrong. The best growth
point on 6TD3 sits beside two arginines and the tool was recommending a basic
group there. With the cutoff it asks for an acceptor, which is what that
environment wants.

**A hypothesis this disproved.** The earlier finding that cavity *fields*
survive docking error far better than cavity *shape* looked like it might be
an artefact of exactly this: a diffuse long-range background does not move
when you nudge the partner, so a field dominated by it would be robust for
trivial reasons. Re-running the sensitivity with the cutoff says no. At 2 Å
of partner error the field correlation goes from 0.915 to 0.907, and at 5 Å
from 0.774 to 0.762. The fields really are more robust than the shape.

Nothing in the gate depends on the electrostatic term, so its thresholds are
unaffected.

### Do the growth points point where chemists actually grew?

Step 6 ranks places on the ligand worth growing chemistry from, driven by how
exposed each atom is once the ligand sits in its pocket. That is a claim, and
the CDK12 series tests it: eight ligands, one purine core, bound in the same
site of the same protein across eight depositions. They are a real
medicinal-chemistry campaign, so wherever they differ from each other is by
definition where chemistry was grown.

`python3 sidecar/run_growth.py`

Superposing the proteins and matching ligand atoms within 1 A splits RC8 into
a core conserved across the whole series -- the purine, C4/C5/C9/N7/C1'/C4'
all 8 of 8 -- and a periphery that varies. Exposure predicts which is which,
with every member of the series used in turn as the reference:

| | median rho | negative |
|---|---|---|
| exposure vs conservation | **-0.50** | 8/8 |
| periphery vs conservation | -0.28 | 8/8 |
| **exposure, periphery removed** | **-0.43** | 8/8, 6 significant |

The third row is the one that matters. A molecule's outer atoms are both more
exposed and more synthetically accessible, so exposure could easily be a proxy
for "not in the middle of the scaffold" and carry nothing the 2D structure does
not already give you. If that were the whole story the partial correlation
would collapse. It does not -- it is stronger than peripherality on its own.

Two honest limits. The eight rows are one series looked at eight ways, not
eight independent tests, so read the consistency rather than a combined
p-value. And this says exposure predicts where a campaign *went*, which is not
the same as saying it predicts where a campaign *should* go: chemists also
grow where the synthesis is easy.

### The ranking is a different claim, and it does not hold

Everything above tests the exposure **gate** -- which atoms get onto the list.
The tool then **orders** that list by a composite of room, direction and
enclosure, and it is easy to let the validated claim lend its credibility to
the unvalidated one. Tested the same way:

| | median rho |
|---|---|
| exposure gate vs conservation | **−0.50** (points at what changed) |
| shipped ranking vs conservation | **+0.37** (points at what did not) |

The ranking puts highest the atoms that campaign left alone. The reason looks
structural rather than like a weighting that needs tuning: conservation
correlates +0.23 with distance to the partner, so the atoms the chemists did
change sit *closer* to DDB1 -- in space a rigid cavity does not show as free.
Real elaboration remodels an interface; this measures the room already in one.

So the chemistry screen says both things, and says which is which. Tuning the
weights until the number came out the other way on a sample of one campaign
would be fitting noise and calling it validation.

**How arbitrary is the order?** The cone that measures "room in front" has
three numbers in it chosen by reasoning rather than measurement. Across six
systems, varying the half-angle between 30 and 45 degrees and the reach
between 5 and 8 angstroms moves the order by Kendall tau 0.81 to 0.96, leaves
the top-ranked atom unchanged in 4 to 5 of 6, and changes the top three by
about one member. The off-axis weighting has no effect on the order at all.

So the best one or two points are robust to those choices and the precision
below them is not. Read the top of the list; ignore the rest of the ordering.

### The gate, and what it costs to calibrate one honestly

The assessment screen answers "is this interface worth doing chemistry into"
with no chemistry done. The temptation is to calibrate it on the 14 solved
glue systems, and that produces a gate that passes whatever it is shown:
every example is a positive, so no threshold drawn from them can fail
anything. It would feel informative and be vacuous.

Negatives were already available. A decoy pose of the partner is a plausible
complex in the wrong place, and it has a real interface with a real cavity at
the ligand site -- exactly the thing the gate has to reject. The hardest ones
are those that bury the presented face best, because those are what a docking
run would rank highest and what the gate would actually be asked about. So
the calibration is 14 crystal poses against the 15 best decoys of each, 210
negatives in all.

`python3 sidecar/run_gate.py`

How well each metric separates the right pose from the wrong ones:

| metric | AUC | in the verdict |
|---|---|---|
| presented face buried | 0.995 | yes, 23% of the weight |
| lined by both proteins | 0.972 | yes, 22% |
| cavity covers the ligand | 0.825 | yes, 15% |
| cavity volume | 0.792 | yes, 14% |
| enclosure quality | 0.788 | yes, 14% |
| mean enclosure | 0.753 | yes, 12% |

All six clear the AUC 0.70 bar, but they are not equal and counting them
equally would let four mediocre metrics outvote two excellent ones. Each is
weighted by how much better than chance it is.

Two details that took a correction each. **Thresholds are the Youden point,
not the weakest validated system.** Taking the weakest as the bar lets one
unusual system drag it down until the hardest decoys clear it, which is what
happened first time: the best wrong poses were scoring borderline. **Where
the two distributions overlap, the check says so** rather than rounding to a
verdict on which side of an arbitrary line a value fell.

What comes out, run back over its own calibration data:

| | pass | borderline | fail |
|---|---|---|---|
| 14 crystal poses | 11 | 3 | **0** |
| 210 wrong poses | **0** | 17 | 193 |

No false passes on 210 adversarial negatives and no false fails on 14
positives. The three borderlines are exactly the three groove-bound ligands,
which is the honest cost of this gate: it would make you hesitate over three
real, published glues. Their sites genuinely are shallower and less enclosed
than the rest, and a gate that waved them through would be waving through the
decoys with them.

### Growing the molecule, and one threshold set from the structures

The chemistry layer is split so the part needing RDKit is as small as it can
be. Scoring a product is numpy and takes coordinates; everything RDKit
touches sits in one function. That boundary is where the tested work stops --
RDKit would not install here either, so the scoring half is tested and the
building half is not.

**What is scored, and what deliberately is not.** Products are judged on
whether they fit, whether they fill the cavity, and whether they reach the
partner. They are not ranked on electrostatic complementarity. The field is a
Coulomb sum over formal charges with no explicit solvent, good enough to say
a region leans one way and nowhere near good enough to order two real
compounds. It is reported beside each product and kept out of the ranking,
the same way cavity volume is measured and kept out of the gate's verdict.

**The clash threshold came from the structures, not from a textbook.** It
started at 2.9 angstroms, which sounds careful and rejected four of the
fourteen real crystal ligands as impossible -- real ligands make hydrogen
bonds at 2.5 to 2.9, and those are contacts. The tightest contact anywhere in
the benchmark is 2.47 A, so the threshold is 2.4 and a test asserts that every
crystal ligand still passes its own test.

| threshold | real ligands rejected |
|---|---|
| 2.4 Å | 0 of 14 |
| 2.5 Å | 1 of 14 |
| 2.9 Å | 4 of 14 |

### What the built-in engine can and cannot do

The engine in `sidecar/glueforge/fftdock.py` exists because LightDock does not
install on Windows. It is Katchalski-Katzir correlation docking: both proteins
on one lattice, the receptor drawn as a scoring shell over a penalised
interior, and one transform pair scoring every translation at once, so the
search is over orientations only.

Two departures from the textbook version, both from what this tool is for.
The ligand is drawn into the receptor grid, so a pose that settles into its
own pocket is penalised while it is generated rather than filtered out
afterwards — which is a bug `dock.py` has to work around on the LightDock
path, because LightDock is given protein only. And the presented face is a
hard geometric filter: a placement is kept only if some partner atom comes
within contact of an exposed ligand atom.

**The scoring works. The search is the limit.** Recognising the right answer
and finding it are different abilities that fail for different reasons, and
a recovery rate blurs them together, so they are measured apart by
`sidecar/run_engine_score.py`: score the crystal placement on the engine's
own grid, score sixty plausible wrong placements the same way, and see where
the crystal one lands.

| entry | crystal | best decoy | recognised | best pose found |
|---|---|---|---|---|
| 6H0F | 72 | 37 | yes | 2.48 Å — hit |
| 6Q0R | 65 | 30 | yes | 8.09 Å |
| 6Q0W | 34 | 27 | yes | 6.08 Å |
| 6Q0V | 26 | 30 | **no** | 7.43 Å |
| 6UD7 | 109 | 26 | yes | 7.68 Å |

6UD7 is the one to read. The grid ranks the crystal pose at four times the
best wrong placement, so the scoring function is not what loses poses —
worth knowing before retuning grid parameters, because the recovery rate on
its own points straight at them. 6Q0V is the exception and no amount of
searching helps there.

**Where it actually stands, measured on the six.** Recovery is 1 of 6 at the
5 Å cutoff. That single number hides two separate failures of similar size,
which `sidecar/compare_rankings.py` pulls apart:

```
best pose available    2.5   8.1   6.1   7.4   7.7  16.1 A   <- the search
best ranking rule      2.6  20.2  20.0  24.9  22.5  34.5 A   <- the choice
```

The search gets within 6–8 Å on four systems and then stops — close, and
not close enough, since the cavity at a glue interface tolerates about 2 Å
of partner error. The ranking then costs another 14 Å of median regret on
top: on five of six systems it picks a pose far worse than one it already
had. Seven ranking rules were compared and the best is barely ahead of the
worst, so this is not a matter of choosing a better formula among the ones
tried.

Both need work and neither alone would be enough. That is a more useful
statement than "17%", and it took three harnesses to be able to make it.

Measured the same way, by rotating the crystal partner a known amount and
taking the best score any translation could then reach:

| off native | 0° | 3° | 5° | 8° | 12° | 20° |
|---|---|---|---|---|---|---|
| best score | +190 | +126 | +136 | +104 | −291 | nothing |

The fit does not decay, it dies. Past about eight degrees no placement near
the right one scores above zero at all, so the search goes elsewhere and
returns something wrong and plausible. That one measurement explains every
failure below, and it is why the presets sample thousands of orientations
and then revisit the best neighbourhoods rather than sampling uniformly.

Two things make a pair hard, and they are independent:

- **A shallow ligand.** A ligand that is not really in a pocket presents a
  broad flat face, and a broad flat face restrains nothing — many wrong
  placements cover it. This was already known from the decoy work (see *Does
  the presented-face restraint actually carry signal?*) and it hurts this
  engine more than LightDock, because the restraint is most of what makes
  this search affordable in the first place.
- **A large partner.** Angular error is set by how finely orientations can be
  sampled, and a large body swings further for the same angle. The
  benchmark's 6500-atom partners are past what the sampling can resolve.

`dock.py` warns on both before spending your time. The measured numbers are
in `cache/recovery/builtin.json`; regenerate with
`python3 sidecar/run_builtin_recovery.py`.

Where LightDock installs, use LightDock. This engine is what you get when it
does not, and it is better than nothing by a margin that depends entirely on
which of the two problems above your pair has.

### What still has to be measured on a machine with LightDock

Everything above is the restraint. How often a restrained *search* finds the
crystal pose is a different number and needs LightDock, which does not
install in the sandbox this was built in. `sidecar/run_recovery.py` runs that
benchmark; until its output exists the three presets show their run times and
say their recovery rates are unknown, which is better than showing a guess.

### What cross-linking MS would add

Less than you would hope for pose accuracy, more than nothing for partner
identification. A 30 A Ca-Ca restraint is about a third the diameter of a small
protein, so it is a proximity statement rather than a geometric constraint.
Starting near the right pose, 5 crosslinks move "within 2 A" from about 40% to
55-70%. In a global search the enrichment is unstable and often indistinguishable
from the shape filter alone.

The restraint you already have is better: the partner must contact the ligand's
presented face, which for CR8 is 6 of 32 heavy atoms. That gives 1.2x to 44x
enrichment depending on the system, median about 7x, on far better sampling than
any crosslink condition.

Two practical notes. Small partner domains are close to uncrosslinkable -
IKZF1's zinc finger has 2 lysines, RBM39's RRM2 has 4, and those are exactly the
substrate classes glues target. And the most valuable use of XL-MS here may not
be restraints at all: differential cross-linking of the glued against the apo
complex answers whether the glue induced a complex, which is a mechanism
question nobody has published for a glue ternary complex.

---

## Getting set up

The landing page does this for you. By hand:

```bash
pip install -r sidecar/requirements.txt    # numpy and scipy, nothing else
pip install lightdock                      # optional: docking
pip install rdkit                          # optional: chemistry
python3 sidecar/fetch.py                   # structures into cache/raw
```

Node is needed only to build and test the page; it has no npm dependencies.

**From a clone, the page has to be built once.** `dist/` is generated, so
it is not in the repository. Starting the server builds it for you if Node
is on PATH, and says so; if Node is not there, the terminal says what to
do before anything opens a browser:

```bash
git clone <this repo> && cd glueforge
python3 glueforge.py     # or double-click GlueForge.bat / .command
```

Node is needed only for that one step, and only on a clone — a release zip
ships the page already built and needs nothing but Python.

The structures in `cache/raw` are not in the repository either — they are
downloadable, so `sidecar/fetch.py` fetches them rather than the history
carrying 14 CIF files. Everything the app needs to run *is* committed: the
five worked examples and the evaluated dashboard data. `python3 glueforge.py
--selftest` passes on a bare clone, and its header says which checks it
skipped for want of structures.

### The worked examples

Five solved ternary complexes ship in `cache/examples`, offered on the
upload screen before anything has been uploaded. Each is the two chains that
matter plus the ligand, in one file, with the anchor chain already
identified — so loading one is a click away from the gate, the cavity view
and the 3D structure.

| entry | ligand | face | buried | cavity | verdict |
|---|---|---|---|---|---|
| 6H0F | Y70 | 72.6 Å² | 83% | 140 Å³ | pass |
| 5FQD | LVY | 74.8 Å² | 82% | 285 Å³ | pass |
| 8BU5 | RMF | 140.7 Å² | 79% | 614 Å³ | pass |
| 6UD7 | EF6 | 168.9 Å² | 68% | 1274 Å³ | pass |
| 6Q0W | EF6 | 343.4 Å² | 36% | 554 Å³ | **borderline** |

Chosen to differ from each other. The benchmark holds eight CDK12 entries
and four DCAF15 ones; shipping them all would be a longer list showing less.
6H0F and 5FQD are the same ligase and chemistry recruiting different
proteins, which is the premise of molecular glue. 8BU5 involves no ligase
pocket at all, so it tests whether the gate measures interfaces rather than
recognising cereblon. 6UD7 has the largest cavity in the set.

**6Q0W is included because it does not pass.** It is the same two proteins
as 6UD7 with a different ligand, and the verdict follows the chemistry
rather than the proteins — indisulam lies in a shallow groove, presenting
343 Å² of broad flat face at 36% buried, and a face like that restrains a
search weakly. A set of five passes would look better and teach less, and
`test_the_example_set_is_not_all_one_thing` keeps it that way.

Every number on every card is produced by `sidecar/make_examples.py` running
the same code the tool runs, and a server test loads each example and checks
the live measurement against what its card claims. Those are different code
paths making the same assertion; if they ever disagree, a card is saying
something the tool contradicts, and the examples are what every later
judgement gets calibrated against.

To rebuild them: `python3 sidecar/make_examples.py` (needs `cache/raw`).

### Which docking engine runs

Docking no longer requires LightDock. If LightDock is on PATH it is used;
if it is not, the built-in engine runs instead and says so.
`dock.py --engine builtin` forces the built-in one even where LightDock
exists.

The default is not a verdict on which is better. The two have never been run
against each other here: `cache/recovery/builtin.json` measures the built-in
engine on six bound-bound cases, and there is no matching file for LightDock
because no machine in this project has had it installed. LightDock is the
default where it exists on the strength of its published benchmarks, not on
anything measured in this repository. Running `sidecar/run_recovery.py
--all-presets` on a machine that has it would produce the comparison; until
then the UI describes the two engines rather than ranking them.

### The built-in docking engine

`sidecar/glueforge/fftdock.py` docks without LightDock, using numpy and
scipy only, so it runs wherever GlueForge does. It is Katchalski-Katzir
correlation docking — both proteins on one lattice, one transform pair
scoring every translation at once — with the presented face as a hard
geometric filter and a local refinement afterwards.

Its measured recovery is in `cache/recovery/builtin.json`, produced by
`sidecar/run_builtin_recovery.py`. Read that before trusting it on anything,
and read `### What the built-in engine can and cannot do` below for where it
falls over. LightDock remains the better engine where it can be installed.

### LightDock does not install on Windows

Not a toolchain problem you can solve by installing a toolchain. LightDock
publishes no Windows wheel for any version, and its own page says plainly:
*"Microsoft Windows is not officially supported, despite many parts of the
protocol might be able to run."* Its `fastdfire` scoring function is compiled
C. Also worth knowing: `freesasa`, which it depends on, has Windows wheels
only up to Python 3.11, so dropping to 3.12 fixes nothing either.

One command is worth trying first, in case a pre-release has changed this:

```
py -m pip install --pre --only-binary=:all: lightdock
```

If that finds nothing, use WSL, where LightDock builds the way its authors
build it:

```
wsl --install -d Ubuntu          # once, then restart
cd /mnt/c/path/to/GlueForge
bash setup-wsl.sh
```

`setup-wsl.sh` installs the toolchain, builds a virtual environment in your
Linux home directory (not on `/mnt/c`, where every scipy import would cross
the filesystem boundary), installs the packages and starts the server. WSL
forwards localhost, so the address it prints opens in your normal Windows
browser, and the project stays where it is.

Docking is the only stage that needs LightDock. Assessing a complex you
already have and growing chemistry into it need numpy, scipy and RDKit, all
of which have Windows builds.

## Evaluating a new pair

The benchmark is not the point of the tool, it is how the tool was checked.
To run it on something of your own:

```bash
python3 sidecar/evaluate.py my_complex.pdb        # mmCIF or PDB, either way
python3 sidecar/evaluate.py 7XYZ --fetch          # or pull it from the PDB
python3 sidecar/evaluate.py my_complex.cif --ligand LVY     # name the glue
python3 sidecar/evaluate.py my_complex.cif --anchor A --partner C
```

It assigns the roles, runs frustration with its gradient, detects the cavity
and its fields, finds the growth points, and prints a report saying whether
the pair follows the benchmark pattern. `--add-to-dashboard` copies the
structure into `cache/raw` so `export_web.py` and `build.js` pick it up and it
appears alongside the others.

**What it needs from you is a structure in which the two proteins are already
positioned relative to each other, with the ligand bound.** That is the one
thing the pipeline cannot work out for itself, and the honest limit of what is
built here:

| what you have | what happens |
|---|---|
| A ternary complex: both proteins, ligand bound | Works today. This is the supported case |
| Anchor + ligand, and the partner, as separate structures | Needs docking first. Not built; the cavity tolerates about 2 A of partner error, which sets the accuracy you need |
| Sequences only | Needs structure prediction, then docking. Neither is here |

A model is as good as a crystal structure to this pipeline, so an AlphaFold or
Boltz ternary prediction works as input - its accuracy just becomes your error
budget.

## Day to day

```bash
python3 sidecar/run_oracle.py              # the frustration hypothesis test
python3 sidecar/run_cavity.py --sensitivity  # cavities and the docking error budget
python3 sidecar/run_restraints.py          # what XL-MS restraints would buy
python3 sidecar/export_web.py              # the frozen cache the dashboard reads
node build.js                              # dist/glueforge.html (served)
node build.js --standalone                 # dist/glueforge-standalone.html
node check.js                              # fast sanity check before a build
npm test                                   # cross-engine + render checks
python3 sidecar/tests/test_engine.py       # engine tests (pytest optional)
python3 sidecar/tests/test_face.py         # presented face + the H convention
python3 sidecar/tests/test_server.py       # the HTTP layer, on a real socket
```

There are two builds because they have opposite jobs. **`glueforge.html`** is
what `glueforge.py` serves: it fetches the validation set from the server the
first time someone opens it, so starting up does not cost four megabytes
nobody asked for. **`glueforge-standalone.html`** is the export, with the data
inlined — no server, no network, no install. Both are checked, because
checking one lets the other drift into the wrong shape unnoticed.

Useful flags:

| command | what it does |
|---|---|
| `run_oracle.py --gamma path/AWSEM_2015.json` | use the published tables instead of the stand-in |
| `run_oracle.py 6TD3 6Q0W` | named entries only |
| `run_oracle.py --decoys 2000` | more decoys per contact (default 800) |
| `run_oracle.py --proximal 8.0` | widen what counts as glue-adjacent (default 6 Å) |
| `fetch.py --list` | show the benchmark set without downloading |
| `build.js --standalone` | inline the cache for an export that needs no server |
| `build.js --out=PATH` | write the build somewhere else |
| `glueforge.py --port NNNN` | pick a port (`--port=NNNN` too); it walks forward from 8731 if that one is taken |
| `glueforge.py --no-browser` | start the server without opening a page |

---

## Potentials

The frustration index needs residue-pair contact energies. Two sources:

- **`awsem-2015` — the default.** The published tables ship in
  `sidecar/data/AWSEM_2015.json`; see `sidecar/data/PROVENANCE.md` for the
  file format and what is and is not implemented. Contact and burial terms;
  the electrostatic term is not implemented.
- **`builtin-v1` — a transparent stand-in**, reachable with
  `--gamma builtin`. Assembled from Kyte–Doolittle hydropathy, formal
  sidechain charge and hydrogen-bond capability, with no burial term. It
  exists so the pipeline runs without the published file and so the
  dependence on the potential can be measured rather than assumed.

Everything else follows the published AWSEM form: `r_min` 4.5 Å, `r_max`
6.5 Å, second shell 6.5–9.5 Å, `eta` 5.0, `eta_sigma` 7.0, `rho_0` 2.6,
`k_contact` 4.184 kJ/mol, CB distances with CA standing in for glycine.

Sign convention, which is easy to get backwards:

```
F = (mean(E_decoy) - E_native) / std(E_decoy)
```

Large positive F is minimally frustrated. Negative F is frustrated. The
conventional cut-offs are +0.78 and −1.0.

---

## How it is put together

```
sidecar/              Python. Everything heavy.
  glueforge/
    cif.py            mmCIF and PDB readers — streaming, stdlib only, no gemmi
    residues.py       published scales and chemical facts, one place
    potential.py      AWSEM form; published tables or the stand-in
    prng.py           mulberry32, bit-identical to the JS side
    frustration.py    the reference engine
    interface.py      anchor/partner assignment and glue labelling
    metrics.py        comparisons, burial stratification, enrichment
    cavity.py         composite cavity detection and pharmacophore fields
    geometry.py       SASA, exit vectors, reach, secondary structure
  evaluate.py         run the whole analysis on one new pair
  fetch.py            RCSB ingest
  run_oracle.py       the frustration hypothesis test
  run_cavity.py       cavities, fields, and the partner-error budget
  run_restraints.py   what XL-MS restraints would and would not buy
  export_web.py       writes cache/web, including the cross-validation fixture
  tests/

src/                  The dashboard. Edit here, not in dist/.
  index.html          shell with three injection markers
  manifest.json       module load order — order matters
  js/
    00-core.js        namespace, formatting, statistics, colour scale
    03-frustration.js the same engine, in the browser, running live
    04-viewer.js      canvas structure viewer
    05-pages.js       benchmark, system, method
    06-cavity.js      the cavity page: fields and the error budget
    50-startup.js     hash routing and boot
  styles/app.css

cache/
  raw/                mmCIF inputs (gitignored)
  web/                the frozen cache the dashboard reads
  oracle/             run_oracle.py output

tests/                Node: cross-engine and render checks
dist/                 generated; not in git
```

### Three tiers

The sidecar does the chemistry and writes everything into a versioned cache.
The browser reads the cache. That split is what lets one build be both a
working tool and a portable demo.

The exception is the frustration engine, which recomputes **in the page** —
strip the glue, change the decoy count, switch flavour, and the map redraws.
That is possible because the frustration index is a statistical potential over
a contact map rather than a simulation: a few hundred lines of arithmetic, not
a dependency.

### The engines are checked against each other

`tests/crossval.test.js` runs the browser engine against a fixture of the
Python engine's outputs. It is exact rather than statistical, because both
draw decoys from the same mulberry32 stream in the same order. Current
agreement is 5 × 10⁻¹⁰, which is the fixture's 9-decimal rounding and not a
real difference.

That check is the reason the dashboard can claim the live engine is validated
rather than merely asserting it.

---

## What this does not do

**The oracle arm only.** Interface geometry comes from the structure that
already contains the answer. That tests whether the hypothesis is true at all,
not whether a pipeline could find it blind. The prospective arm — blind A–B
docking, partner recovery, pocket recovery — is not built.

**Frustration never touches the ligand.** The index is defined over
amino-acid pairs; there is no small-molecule term in the Hamiltonian. The
ligand only labels which contacts count as proximal. Any phrase like "this
compound relieves that contact" is an inference through geometry.

**All four systems share DDB1.** CDK12–cyclin K recruits it directly; CRBN and
DCAF15 are both DDB1 substrate receptors. A signal common to all of them may
be about DDB1-family interfaces rather than glue interfaces generally. Adding
a 14-3-3 stabiliser system would break that confound and is the single highest
-value addition to the benchmark.

**Other cut-offs are still untested.** The proximity cut-off has been swept
from 4 to 10 Å and the threshold-free correlation agrees with it. The contact
cut-off (9.5 Å) and the decoy count have not been varied.

**No chemistry.** Enumeration, scoring and the enrichment metric are
specified in the proposal and not implemented here.

**No docking.** The cavities are computed from crystal ternary complexes, so
the partner is placed perfectly. The error budget says what happens when it is
not, but the search that would produce a real placement needs LightDock or
Boltz-2 and runs on your machine.

**Frustration SAR does not work on crystal series.** Across the eight CDK12
glues, between-ligand variation in per-residue frustration is 0.151 against a
within-structure noise floor of 0.114, and the residues that vary most are the
ones furthest from the ligand. That variation is crystallography, not chemistry.
The in-silico version - one structure, swapped ligands, minimised - is the
clean test and needs OpenMM.

---

## Where to go next

1. **Re-run with `AWSEM_2015.json`.** One flag. Until then every frustration
   number is about the stand-in. The cavity results do not depend on it, which
   is one more reason that side is on firmer ground.
2. **Add a 14-3-3 system** to break the DDB1 confound.
3. **Build step 4's scoring half** — a ligand-aware interface filter, validated
   by generating decoy poses and checking it ranks the crystal pose top. That
   is what lets a new protein pair be brought in rather than only re-examined.
4. **Build the prospective arm** — the expensive half.
