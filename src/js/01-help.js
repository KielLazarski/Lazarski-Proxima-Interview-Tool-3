/* GlueForge — the glossary and its popover.
 *
 * Every number on a page can be asked what it means, at three depths. The
 * three are written to be genuinely different registers rather than the same
 * sentence at three lengths:
 *
 *   eli5     — no jargon at all, an analogy if one helps
 *   concise  — what a colleague would say in a corridor
 *   detail   — the formula, the units, the caveat that bites
 *
 * Terms are keyed by a short slug and attached to an element with
 * GF.help.attach(node, 'slug').
 */

(function (GF) {
  'use strict';

  const TERMS = {

    /* ---- the index itself --------------------------------------------- */

    frustration: {
      title: 'Frustration index (F)',
      eli5:
        'Take two amino acids that are touching. Now imagine swapping them ' +
        'for a random pair. If most random pairs would do the job just as ' +
        'well or better, that contact is "frustrated" — the protein is ' +
        'putting up with it rather than benefiting from it. If the real pair ' +
        'beats almost every alternative, it is not frustrated.',
      concise:
        'A z-score. For each contact we compute its energy, then the energy ' +
        'of several hundred decoy pairings at the same site, and ask how many ' +
        'standard deviations better the real one is. Positive means the ' +
        'native pair beats its alternatives; negative means it does not.',
      detail:
        'F = (mean(E_decoy) − E_native) / sd(E_decoy), with energies in the ' +
        'AWSEM contact form over CB–CB pairs (CA for glycine). Conventional ' +
        'cut-offs: F > +0.78 minimally frustrated, F < −1.0 highly ' +
        'frustrated. Note the sign — a lower native energy is more ' +
        'favourable, so it appears as a larger positive F. The index is ' +
        'defined over amino-acid pairs only; there is no small-molecule term, ' +
        'so the ligand never enters this number.'
    },

    interface_contacts: {
      title: 'Interface contacts',
      eli5:
        'The places where the two proteins actually touch each other. Not ' +
        'the whole protein — just the handshake.',
      concise:
        'Residue pairs within 9.5 Å where one residue belongs to the anchor ' +
        'protein and the other to the partner. Contacts inside a single ' +
        'protein are computed too, but excluded from this count.',
      detail:
        'A contact is a CB–CB pair under the 9.5 Å cut-off, which is where ' +
        'the mediated well of the potential goes to zero. Pairs fewer than 2 ' +
        'apart in sequence are dropped so backbone neighbours do not fill the ' +
        'list; that rule only applies within a chain, since sequence ' +
        'separation is undefined across one.'
    },

    proximal_distal: {
      title: 'Glue-proximal vs distal',
      eli5:
        'We split the handshake into two piles: the contacts right next to ' +
        'where the drug sits, and all the rest. Then we compare the two ' +
        'piles. That comparison is the whole experiment.',
      concise:
        'Interface contacts within 6 Å of any ligand heavy atom are ' +
        '"proximal"; the remaining interface contacts are "distal". A ' +
        'contact counts as proximal if either of its two residues reaches ' +
        'the ligand.',
      detail:
        'The cut-off is one solvation shell beyond direct contact and is ' +
        'adjustable with --proximal. No sensitivity sweep over it has been ' +
        'run, which is a real gap: if the sign flips between 4 Å and 10 Å ' +
        'the result is about the cut-off rather than about the chemistry. ' +
        'Where the interface is small and glue-dominated, the distal group ' +
        'gets thin and the comparison loses power — 6H0F has only 27 distal ' +
        'contacts.'
    },

    difference: {
      title: 'Difference (Δ)',
      eli5:
        'Mean frustration near the drug, minus mean frustration everywhere ' +
        'else on the interface. A negative number means the contacts by the ' +
        'drug are the more frustrated ones — which is what we predicted.',
      concise:
        'mean(F_proximal) − mean(F_distal). Negative supports the ' +
        'hypothesis that a glue sits where the interface is frustrated.',
      detail:
        'An unstandardised difference in a z-score, so it is already in ' +
        'standard-deviation-ish units of the decoy distribution, but it is ' +
        'not an effect size — compare it against the rank-biserial, which ' +
        'is. The magnitudes here run roughly −0.6 to +0.15.'
    },

    burial_controlled: {
      title: 'Burial-controlled difference',
      eli5:
        'A glue sits in the middle of the handshake, and the middle of any ' +
        'handshake is different from its edges for reasons that have nothing ' +
        'to do with our idea. So we redo the comparison inside groups of ' +
        'contacts that are buried to a similar degree, and see if the ' +
        'difference survives.',
      concise:
        'The same proximal-vs-distal difference, computed within quartiles ' +
        'of mean contact density and pooled by band size. Controls for a ' +
        'glue occupying the interface core.',
      detail:
        'Covariate is 0.5·(rho_i + rho_j), the AWSEM local density averaged ' +
        'over the contact. Four quantile bands, difference computed per band, ' +
        'pooled weighted by band n. Across the current benchmark this shrinks ' +
        'the mean effect by roughly 20% without removing it, so the raw ' +
        'result is not purely a core-versus-rim artefact — though a four-band ' +
        'stratification is a coarse control, not a regression.'
    },

    rank_biserial: {
      title: 'Rank-biserial correlation',
      eli5:
        'If you picked one contact from each pile at random, how often would ' +
        'the one near the drug be the more frustrated of the two? This number ' +
        'is that tendency, from −1 to +1. Zero means a coin flip.',
      concise:
        'A non-parametric effect size derived from the Mann-Whitney U. ' +
        'Unlike the raw difference it is bounded and insensitive to the ' +
        'shape of the distribution, which matters because frustration is not ' +
        'normally distributed.',
      detail:
        'r = 1 − 2U/(n_a·n_b). Here positive values mean the proximal group ' +
        'ranks lower in frustration index, i.e. more frustrated. Values in ' +
        'this benchmark run +0.19 to +0.46 where the effect is present, which ' +
        'is small-to-moderate.'
    },

    p_permutation: {
      title: 'p (permutation)',
      eli5:
        'How often you would see a gap this big purely by chance, if the two ' +
        'piles were really the same. Small means "probably not luck".',
      concise:
        'Two-sided permutation test on the difference in means, 20,000 ' +
        'resamples. No distributional assumption, which is the point — the ' +
        'frustration index is skewed.',
      detail:
        'These p-values are uncorrected across fourteen structures. They are ' +
        'also not independent: eight of the fourteen are the same interface ' +
        'with different ligands. Treat any single row as descriptive and read ' +
        'the across-system sign test on the benchmark page instead.'
    },

    density: {
      title: 'Local density (rho)',
      eli5:
        'How crowded it is around a residue. Buried in the core of a protein ' +
        'is crowded; sticking out into water is not.',
      concise:
        'The AWSEM local density: a smooth count of neighbours in the ' +
        'direct-contact window. It decides whether a mediated contact is ' +
        'treated as water-mediated or protein-mediated.',
      detail:
        'rho_i = sum over j of theta(r_ij), with theta the 4.5–6.5 Å double-' +
        'tanh window and pairs closer than 2 in sequence excluded. The ' +
        'potential switches on rho_0 = 2.6; the observed mean across these ' +
        'interfaces is about 2.6, which is a useful sign the implementation ' +
        'is behaving.'
    },

    anchor_share: {
      title: 'Anchor share',
      eli5:
        'The drug is mostly buried in one of the two proteins and only pokes ' +
        'a small face out at the other. This is how much goes into the first ' +
        'one. It is also how the software worked out which protein is which, ' +
        'with nothing configured by hand.',
      concise:
        'Fraction of protein heavy-atom contacts within 4.5 Å of the ligand ' +
        'that belong to the anchor chain. The chain burying more of the ' +
        'ligand is taken as the anchor; the next one as the partner.',
      detail:
        'Running 0.68–0.84 across this benchmark, which independently ' +
        'reproduces the asymmetry reported for dWIZ (CRBN–ZF7 414 Å² of ' +
        'buried surface against 152 Å² for the glue contact). It is the ' +
        'structural basis for the claim that glue design is mostly pair ' +
        'selection: a 400–500 Da molecule has nowhere near enough surface to ' +
        'build the interface itself.'
    },

    /* ---- controls ------------------------------------------------------ */

    decoys: {
      title: 'Decoys per contact',
      eli5:
        'How many random stand-in pairs each real contact is measured ' +
        'against. More stand-ins give a steadier score and take longer. 800 ' +
        'is plenty; 200 starts to look noisy.',
      concise:
        'Sample size for the decoy energy distribution behind each z-score. ' +
        'The mean and standard deviation of that sample are the denominator ' +
        'of F, so too few decoys inflates the variance of the index itself.',
      detail:
        'Cost is linear in decoys × contacts. The decoys are drawn from a ' +
        'seeded mulberry32 stream, so a given setting is reproducible to the ' +
        'bit — the same button press always gives the same answer, and the ' +
        'Python and JavaScript engines give the same answer as each other.'
    },

    flavour: {
      title: 'Configurational vs mutational',
      eli5:
        'Two ways of building the random stand-ins. "Mutational" only swaps ' +
        'which amino acids are there and keeps everything where it is. ' +
        '"Configurational" also shuffles how far apart and how buried they ' +
        'are. Configurational is the harder test to pass.',
      concise:
        'Decoy generation mode. Mutational resamples residue identity alone, ' +
        'holding distance and burial at their native values. Configurational ' +
        'resamples identity, distance and local density, each drawn from the ' +
        'structure\'s own distributions.',
      detail:
        'Mutational answers "is this the right pair of residues for this ' +
        'spot". Configurational answers "is this a good contact at all, ' +
        'compared with the contacts this protein makes generally". The ' +
        'benchmark reports configurational because the question is about ' +
        'interface quality rather than sequence optimality; both are cached ' +
        'and the sidecar writes either.'
    },

    recompute: {
      title: 'Recompute in page',
      eli5:
        'Runs the whole frustration calculation again, right here in the ' +
        'browser, instead of reading the answer that was computed earlier. ' +
        'It is the one part of this pipeline fast enough to do that.',
      concise:
        'Re-runs the JavaScript engine on the cached coordinates with the ' +
        'settings you chose, then reports how long it took and how far it ' +
        'landed from the cached run.',
      detail:
        'Everything else on these pages is read from a frozen cache written ' +
        'by the Python sidecar. Frustration recomputes live because it is a ' +
        'statistical potential over a contact map — arithmetic and table ' +
        'lookups — rather than a simulation. Changing decoys or flavour and ' +
        'getting a different number is the honest demonstration that the ' +
        'engine is really running and not replaying.'
    },

    histogram: {
      title: 'Frustration against distance to the ligand',
      eli5:
        'Every dot is one point of contact between the two proteins. How ' +
        'far left it sits is how frustrated that contact is. How low it sits ' +
        'is how close it is to the drug. The shaded strip at the bottom is ' +
        'the group right next to the drug. If the idea behind this whole ' +
        'project is right, that strip should lean further left than the rest.',
      concise:
        'Frustration index on the horizontal axis against minimum distance ' +
        'from the contact to any ligand heavy atom on the vertical. The ' +
        'dashed horizontal line is the proximity cut-off that defines the ' +
        'two groups; vertical rules mark each group\'s mean.',
      detail:
        'An earlier version plotted the two groups as jittered rows, which ' +
        'put nothing on the vertical axis and invited it to be read as data. ' +
        'Distance is the variable that actually defines the split, so ' +
        'plotting against it shows whether the effect is a genuine gradient ' +
        'or an artefact of where the cut-off was drawn \u2014 the open ' +
        'question flagged on the method page. A contact\'s distance is the ' +
        'smaller of its two residues\' distances to the ligand.'
    },

    gradient: {
      title: 'Gradient toward the ligand',
      eli5:
        'Instead of splitting the contacts into two piles and comparing ' +
        'them, this just asks: as you walk toward the drug, does frustration ' +
        'go up? A positive number means yes. It avoids having to pick a ' +
        'dividing line at all.',
      concise:
        'Spearman rank correlation between a contact\'s frustration index ' +
        'and its distance to the ligand, over interface contacts. Positive ' +
        'means frustration rises as the ligand is approached. The sweep ' +
        'recomputes the group difference at 4, 5, 6, 8 and 10 Å.',
      detail:
        'This exists because the split comparison cannot, on its own, tell a ' +
        'real gradient from a cut-off that happened to land somewhere ' +
        'flattering. Across this benchmark rho is positive in 13 of 14 ' +
        'structures, mean +0.23, and the sweep keeps its sign in almost all ' +
        'of them, so the 6 Å choice is not driving the result. The two CRBN ' +
        'systems stay flat at every cut-off, which makes that null cleaner ' +
        'than a single threshold comparison could.'
    },

    cavity_volume: {
      title: 'Cavity volume',
      eli5:
        'How much empty space there is at the join between the two proteins ' +
        'for a molecule to fill. Too small and nothing fits; very large and ' +
        'it is probably not one pocket but several run together.',
      concise:
        'Volume of the enclosed free space at the glue site, found by ' +
        'rasterising both proteins and keeping connected points that are ' +
        'free, enclosed and far enough from protein to hold a probe.',
      detail:
        'LIGSITE-style detection on a 0.8 A grid with a 1.4 A probe, plus a ' +
        'morphological closing step to bridge the thin necks that otherwise ' +
        'fragment a shallow groove into dozens of pieces. Validated by ' +
        'checking the real ligand lands inside: 11 of 14 systems recover ' +
        '100% of its heavy atoms.'
    },

    enclosure_quality: {
      title: 'Enclosure quality',
      eli5:
        'Whether the site is a proper pocket with walls around it, or just a ' +
        'dent in the surface. Dents are much harder to design for.',
      concise:
        'Share of cavity points enclosed along all three axes. Near 1 is a ' +
        'buried pocket; near 0 is a surface groove.',
      detail:
        'This matters more than it looks. Pocket-conditioned generative ' +
        'models are trained overwhelmingly on enclosed druggable pockets, so ' +
        'a groove is out of distribution for them. Across this benchmark the ' +
        'DCAF15 Faust form scores 0.05 to 0.25 and is the one site where the ' +
        'detector does not recover the ligand fully \u2014 which matches the ' +
        'published description of it as a shallow, non-conserved pocket.'
    },

    cavity_split: {
      title: 'Anchor-lined share',
      eli5:
        'Most of the pocket is made by the protein the drug is already stuck ' +
        'to, and you have a crystal structure of that. Only the smaller rest ' +
        'depends on guessing where the second protein goes.',
      concise:
        'Fraction of cavity points whose lining comes mostly from the anchor ' +
        'rather than the partner. It runs 61 to 76% across every system here.',
      detail:
        'This is why the design problem is better conditioned than generic ' +
        'PPI design. Under a 2 A partner displacement the anchor-lined region ' +
        'holds up at Jaccard 0.87 to 0.97 while the partner-lined region ' +
        'falls to 0.49 to 0.78. The practical recipe: hard constraints from ' +
        'the anchor side, soft constraints from the partner side.'
    },

    cavity_field: {
      title: 'Pharmacophore fields',
      eli5:
        'Maps of what the pocket wants: where a positive or negative charge ' +
        'would be happy, where it can make hydrogen bonds, and where it is ' +
        'greasy. These are what you design a molecule to match.',
      concise:
        'Per-point electrostatic potential, hydrogen-bond donor and acceptor ' +
        'proximity, hydrophobic contact, and how much of each point\'s lining ' +
        'comes from the anchor. A pocket-conditioned generative model is ' +
        'conditioned on maps like these.',
      detail:
        'Electrostatics are formal charges on ionisable groups with a ' +
        'distance-dependent dielectric (eps = 4r). That gets the acid-or-base ' +
        'character of a cavity right and nothing finer; anything quantitative ' +
        'needs real partial charges and a Poisson-Boltzmann solver. The ' +
        'fields are far more robust to partner displacement than the shape ' +
        'is, which argues for weighting them above shape in a design loop.'
    },

    error_budget: {
      title: 'Docking error budget',
      eli5:
        'If you guess the second protein\'s position slightly wrong, how ' +
        'wrong does the pocket become? The answer sets how good your docking ' +
        'has to be before designing into it is worth doing.',
      concise:
        'The partner is displaced rigidly by a given amount and the cavity ' +
        'recomputed. Jaccard is shape agreement with the undisplaced cavity; ' +
        'the r columns are field correlations over the overlap.',
      detail:
        'Shape is usable to about 2 A and falls apart by 3 A. The fields hold ' +
        'up far better \u2014 electrostatic correlation stays above 0.78 even ' +
        'at 5 A \u2014 because they are dominated by the anchor side and by ' +
        'long-range terms. Everything here is rigid-body with no induced fit, ' +
        'so a real partner relaxing into a displaced position should degrade ' +
        'more gently than this.'
    },

    style: {
      title: 'Drawing style',
      eli5:
        'Four ways of drawing the same two proteins. Contacts is the one to ' +
        'read the science from. The ribbons are for getting your bearings in ' +
        'the fold first. Exit vectors shows where chemistry could grow.',
      concise:
        'Contacts draws a thin backbone with contact chords. Ribbon draws a ' +
        'smoothed CA spline, wider through helices and strands. Side chains ' +
        'adds sticks on interface residues only. Exit vectors overlays the ' +
        'ligand growth points.',
      detail:
        'Secondary structure is estimated from CA i→i+3 and i→i+4 distances, ' +
        'not from DSSP, because the cache carries no backbone N, C or O. It ' +
        'is good enough to draw a fold that reads correctly and should not ' +
        'be quoted as a secondary-structure assignment. Side chains are ' +
        'exported for interface residues only, to keep the build portable.'
    },

    presented_face: {
      title: 'Presented face',
      eli5:
        'Once the drug is sitting in its pocket in the first protein, part ' +
        'of it still sticks out into the water. That sticking-out part is ' +
        'the only bit the second protein can ever touch — and the only place ' +
        'you can grow new chemistry from.',
      concise:
        'Solvent-accessible surface of the ligand computed against the ' +
        'anchor chain alone, not the ternary complex. In the ternary the ' +
        'glue is sandwiched between both proteins and is buried by ' +
        'definition, which would report no growth points at all.',
      detail:
        'Shrake-Rupley with a 1.4 Å probe over 192 golden-spiral points per ' +
        'atom, Bondi radii, heavy atoms only. "Buried by the partner" is ' +
        '1 − SASA_ternary / SASA_binary. Hydrogens are discarded whether or ' +
        'not the structure has them: 9 of these 14 depositions model them and ' +
        '5 do not, and a ligand’s own hydrogens shield its own heavy ' +
        'atoms enough to change its apparent exposed area by 1.7× to 4.3×. ' +
        'Keeping them would make absolute areas incomparable between ' +
        'structures, so they go. The cost is a little physical accuracy on ' +
        'the structures that have them; prefer buried fraction when ' +
        'comparing systems, since it is stable across both conventions to ' +
        'about ±0.05 where absolute area is not. On this benchmark the ' +
        'presented faces run from about 73 Å² to 343 Å², and they are small ' +
        'against protein-protein interfaces of several hundred to a few ' +
        'thousand. That is the structural reason glue design is mostly pair ' +
        'selection.'
    },

    addressable: {
      title: 'Addressable frustration',
      eli5:
        'How much frustration a substituent growing from this atom could ' +
        'plausibly reach. Bigger is a more promising place to grow.',
      concise:
        'Sum over reachable frustrated contacts of |F|, discounted by how ' +
        'far each sits and how far off the exit direction it lies. It ranks ' +
        'growth points against each other; the absolute value has no units ' +
        'worth quoting.',
      detail:
        'Weight is (1 − d/d_max) × normalised angular agreement above the ' +
        'cut-off, with d_max 12 Å and the cone half-angle 75°. Twelve ' +
        'angstroms is roughly what three to five rotatable bonds span. Only ' +
        'contacts with F below −0.3 are counted. This is a feasibility and ' +
        'prioritisation score, not a prediction that any particular molecule ' +
        'will work.'
    },

    exit_vectors: {
      title: 'Exit vectors',
      eli5:
        'Places on the drug you could attach something new, and which way ' +
        'that something would point. The surprise is that the most exposed ' +
        'atom is often not the best one — plenty of them point away from the ' +
        'other protein entirely.',
      concise:
        'Exposed ligand heavy atoms with their local outward normals, taken ' +
        'as the vector from the mean of an atom\'s bonded neighbours to the ' +
        'atom itself. That is where a substituent would project.',
      detail:
        'Atoms below 8% relative exposure are dropped. Bonds are inferred ' +
        'geometrically at 0.4–1.95 Å rather than read from a connectivity ' +
        'table, so an unusual ligand may need checking. Turning a chosen ' +
        'vector into an actual molecule needs RDKit and is not in this ' +
        'build; choosing which vector is worth the effort does not, which is ' +
        'why this part exists.'
    },

    viewer: {
      title: 'The structure view',
      eli5:
        'A simple line drawing of the two proteins. Each one gets its own ' +
        'colour. Dots mark the residues that touch the other protein, and ' +
        'the orange blob is the drug. The coloured lines are the contacts, ' +
        'red where frustrated and blue where not.',
      concise:
        'Orthographic backbone trace, one hue per chain. Filled nodes are ' +
        'interface residues. Chords are contacts coloured on the diverging ' +
        'frustration scale; the glue-proximal ones are drawn heavier.',
      detail:
        'Deliberately a flat line drawing rather than a lit 3D render: the ' +
        'thing being read is which contacts are frustrated and where they ' +
        'sit relative to the ligand, and ribbons get in the way of that. The ' +
        'ligand is drawn for reference and never enters any calculation.'
    }
  };

  /* ---- popover -------------------------------------------------------- */

  let open = null;

  function close() {
    if (open && open.parentNode) open.parentNode.removeChild(open);
    open = null;
  }

  function build(slug) {
    const t = TERMS[slug];
    if (!t) return null;

    const pop = GF.el('div', { class: 'pop', role: 'dialog' });
    pop.appendChild(GF.el('div', { class: 'pop-title', text: t.title }));

    const tabs = GF.el('div', { class: 'pop-tabs' });
    const bodyEl = GF.el('div', { class: 'pop-body' });

    const levels = [
      ['eli5', 'Plain'],
      ['concise', 'Concise'],
      ['detail', 'Detail']
    ];

    function show(key, btn) {
      bodyEl.textContent = t[key];
      Array.from(tabs.children).forEach(function (c) { c.classList.remove('on'); });
      btn.classList.add('on');
      try { localStorage.setItem('gf.help.level', key); } catch (e) { /* private mode */ }
    }

    let preferred = 'concise';
    try {
      const saved = localStorage.getItem('gf.help.level');
      if (saved && TERMS[slug][saved]) preferred = saved;
    } catch (e) { /* ignore */ }

    levels.forEach(function (lv) {
      const b = GF.el('button', { class: 'pop-tab', text: lv[1] });
      b.addEventListener('click', function (e) {
        e.stopPropagation();
        show(lv[0], b);
      });
      tabs.appendChild(b);
      if (lv[0] === preferred) setTimeout(function () { show(lv[0], b); }, 0);
    });

    pop.appendChild(tabs);
    pop.appendChild(bodyEl);
    pop.addEventListener('click', function (e) { e.stopPropagation(); });
    return pop;
  }

  /** Add a small "?" affordance to a node that opens the glossary entry. */
  function attach(node, slug, opts) {
    if (!TERMS[slug]) return node;
    const mark = GF.el('button', {
      class: 'helpmark' + (opts && opts.inline ? ' inline' : ''),
      text: '?',
      type: 'button',
      'aria-label': 'What is ' + TERMS[slug].title + '?'
    });

    mark.addEventListener('click', function (e) {
      e.stopPropagation();
      const wasOpen = open && open.dataset && open.dataset.slug === slug;
      close();
      if (wasOpen) return;
      const pop = build(slug);
      if (!pop) return;
      pop.dataset.slug = slug;
      mark.parentNode.appendChild(pop);
      open = pop;
    });

    node.appendChild(mark);
    node.classList.add('has-help');
    return node;
  }

  /** A label plus its help affordance, as one node. */
  function label(text, slug, tag) {
    const n = GF.el(tag || 'span', { class: 'hlabel' }, [
      GF.el('span', { text: text })
    ]);
    return attach(n, slug);
  }

  document.addEventListener('click', close);
  document.addEventListener('keydown', function (e) {
    if (e.key === 'Escape') close();
  });

  GF.help = { TERMS: TERMS, attach: attach, label: label, close: close };
})(globalThis.GF || (globalThis.GF = {}));
