# Onderzoek: HOIJ implementeren in lavaan

**HOIJ = Higher-Order Infinitesimal Jackknife** (tweede-orde IJ, "HOIJ-2").
Doel van de implementatie: **standaardfouten en betrouwbaarheidsintervallen**
voor (functies van) modelparameters van een gefit lavaan-model, zonder
bootstrap-herfits — de bootstrapverdeling wordt benaderd via een tweede-orde
Taylor-expansie van de schatter in de observatiegewichten.

Bronnen in deze repo:
- `kernfuncties` — de rekenkern (casewise loglik, casewise informatie `J`,
  derde-orde tensor `T`, α-kalibratie).
- `HOIJ_simulatiestudie` — volledige toepassing: de HOIJ-2-replicatielus
  (§8, r. 1108–1180), gedeelde multinomiale gewichten, percentiel-CI's, en
  een pre-flight zelftest die alle schaalconventies verifieert (§6, r. 702–778).

---

## 1. Hoe HOIJ werkt (zoals toegepast in `HOIJ_simulatiestudie`)

Eenmalige setup per gefit model (θ̂ = `coef(fit)`, D parameters, N cases):

1. **Scores** `S` (N×D): `lavScores(fit, scaling = TRUE)` — per case de
   gradiëntbijdrage (conventie: −sᵢ/N, geverifieerd in zelftest (a)).
2. **Inverse geobserveerde informatie** `H⁻¹`:
   `lavTech(fit, "inverted.information.observed")`.
3. **Casewise geobserveerde informatie** `J` (N×D×D): `compute_all_J()` —
   numerieke tweede afgeleiden van de casewise loglik (zelftest (b):
   `information.observed = Σᵢ Jᵢ / N`).
4. **Derde-orde tensor** `T` (D×D×D): `compute_T_tensor_grad()` op lavaan's
   analytische gradiënt van F, teruggeschaald naar loglik-schaal via
   `calibrate_alpha()` (zelftest (c)/(d)).
5. **Gewichten**: B multinomiale countvectoren `W` (B×N),
   `Δw = W − 1` — *dezelfde* gewichten als een exacte bootstrap zou
   gebruiken, zodat HOIJ−bootstrap-verschillen pure approximatiefout zijn.

Per gewichtsvector (rij i van Δw):

```
c   = H⁻¹ Sᵀ Δwᵢ                 (invloedsterm)
IJ1:    θ(w) ≈ θ̂ − c                                  (1e orde)
HOIJ-2: d₂ = H⁻¹ J(Δwᵢ) c  −  ½ H⁻¹ T (c ⊗ c)
        θ(w) ≈ θ̂ − c + s·d₂,   s = min(1, κ·‖c‖/‖d₂‖)  (trust-region-demping, κ = 0.5)
```

Daarna wordt de functionaal φ (bv. `ab`, `speed~~speed`, R², omega) op alle
B pseudo-replicaten geëvalueerd: **SE = sd van de replicaatwaarden,
CI = percentielinterval**. HOIJ-replicaten kunnen per constructie niet
"niet convergeren" — het grote praktische voordeel boven bootstrap
(1 fit + goedkope algebra i.p.v. B herfits).

Diagnostiek die de studie meeneemt en die een implementatie moet behouden:
α-`spread` (kalibratieconsistentie, tolerantie 0.1), dempingsfractie/mean s,
fractie niet-toelaatbare replicaten (negatieve varianties), en een
sensitiviteitsvariant die niet-toelaatbare replicaten schrapt.

---

## 2. Wat lavaan al levert vs. wat zelf gebouwd moet worden

**Al beschikbaar via geëxporteerde API** (geen risico):
- `lavScores(fit, scaling = TRUE)` — casewise scores.
- `lavTech(fit, "information.observed")` / `"inverted.information.observed"`.
- `lavTech(fit, "inverted.information.expected")` — voor de Wald-vergelijking.
- `parTable()`-route voor eventuele verificatie-bootstraps.

**Zelf te bouwen (nu prototype, met verbeterpunten):**

| Onderdeel | Nu | Gewenst |
|---|---|---|
| Casewise `Jᵢ` (N×D×D) | numerieke 2e FD op casewise loglik, O(D²) loglik-evaluaties | analytisch via kettingregel: Jᵢ = Δᵀ Hᵢ(μ,Σ) Δ + 2e-afgeleide-term van de implied moments; lavaan's Δ-Jacobiaan (`lav_model_delta`) levert de eerste helft |
| `T`-tensor (D×D×D) | numerieke 2e FD op analytische ∂F/∂θ + α-schaalbrug | idem numeriek is acceptabel, maar dan met Richardson-extrapolatie/`numDeriv`; α wordt overbodig zodra alles op één (loglik-)schaal staat |
| α-kalibratie | empirische mediaan-ratio + spread-check | degraderen tot pure *diagnostiek* (zelftest), niet tot rekenstap |
| `lavaan:::`-internals (`lav_model_x2glist`, `lav_model_implied`, `lav_model_gradient`) | direct aangeroepen | vervangen door `lav_export_estimation(fit)` waar mogelijk; anders versie-pin + CI-tests |

---

## 3. Voorstel publieke API

Eén gebruikersfunctie (werknaam) die de hele §8-kern van de simulatiestudie
generiek maakt:

```r
hoijLavaan(fit,
           functional = NULL,   # NULL = alle vrije parameters;
                                # of character ("a*b", ":=" -achtige expressie)
                                # of functie phi(theta) / lijst daarvan
           B = 1000L,           # aantal gewichtsvectoren
           order = 2L,          # 1 = IJ1, 2 = HOIJ-2
           kappa = 0.5,         # trust-region-demping
           ci = c("perc", "none"),
           level = 0.95,
           admissibility = c("keep", "drop"),  # cf. hoij2 vs hoij2_sens
           seed = NULL)
```

Retourneert een object met: schatting, **SE**, CI-grenzen, en diagnostiek
(α-spread, dempingsfractie, fractie niet-toelaatbaar, timing) + nette
`print()`/`summary()`. Aansluiting bij lavaan-conventies: zelfde interface-stijl
als `bootstrapLavaan()`; functionalen bij voorkeur via lavaan's `:=`-syntaxis
zodat gebruikers geen R-functies hoeven te schrijven (intern vertalen naar
φ(θ) op de vrije-parametervector, zoals `functionals_scalar/vec` in de studie).

Ontwerpkeuzes die uit de studie meegenomen moeten worden:
- **Gedeelde/reproduceerbare gewichten** (seed-argument; optie om `W` te
  exporteren zodat een gebruiker HOIJ tegen een echte bootstrap kan leggen).
- **Geen stille fallbacks**: NA + reden (fallb_Hobs / fallb_alpha / fallb_deriv)
  in plaats van geruisloos degraderen naar IJ1.
- **Vectorisatie van de B-lus** zoals in het script (C = ΔW·S·H⁻¹ als één
  matrixproduct; J-contractie via `ΔW %*% J_2d`), zodat de kosten na de setup
  O(B·D²) blijven.

---

## 4. Dekkingsgaten t.o.v. lavaan's generaliteit (scopebepaling)

De kern gaat impliciet uit van **één groep + complete data + ML zonder
meanstructure** (overal `[[1]]`; μ-fallback `colMeans(X)`; zie ook de
covmatrix-bootstrap-route die op sufficiency onder ML leunt). Voor een
eerste release is dat een verdedigbare scope, mits hard afgedwongen:

- expliciete checks op `estimator == "ML"`, 1 groep, geen missing, geen
  ordinale variabelen, geen (on)gelijkheidsrestricties — met informatieve
  foutmeldingen;
- meanstructure netjes ondersteunen (θ bevat dan ook intercepten; de
  casewise loglik en scores dekken dat al, alleen de μ-fallback moet weg);
- daarna uitbreiden: multi-group (blok-structuur in S, J, T), FIML/missing
  (casewise loglik per missing-patroon), MLR/weights;
- restricties/`:=`-parameters: `:=`-functionalen zijn al gedekt via φ(θ);
  échte gelijkheids-/ongelijkheidsrestricties vergen de constraint-Jacobiaan
  (raakvlak met restriktor) en zijn fase 2+.

---

## 5. Numeriek & performance

- FD-stapgroottes (`delta = 1e-5`, `h = 1e-4`) zijn nu hard-coded en op één
  model getest; maak ze schaal-bewust (relatief t.o.v. |θₖ|) of gebruik
  `numDeriv::genD`-achtige extrapolatie. De zelftest (a)–(f) uit de studie is
  het juiste vangnet — die hoort als unit-test in het pakket.
- Setupkosten: `J` kost O(D²) casewise-loglik-evaluaties, `T` kost O(D²)
  gradiënt-evaluaties; bij D = 21 is dat ~2×231 evaluaties — prima, maar het
  groeit kwadratisch. Analytische `Jᵢ` (via Δ) is de belangrijkste versnelling
  voor grotere modellen.
- Geheugen: `J` is N×D×D (bij N = 500, D = 21: ~1.8 MB — geen probleem;
  bij D = 100 wel ~40 MB per 500 cases — documenteren, evt. chunken).
- Randgevallen: bijna-singuliere `H` (nu `ginv`-fallback), Heywood-cases in
  replicaten (toelaatbaarheidsdiagnostiek behouden), grote `spread` → fout.

---

## 6. Verpakking, tests en validatie

- **R-pakket** (standalone, `Imports: lavaan`), roxygen-docs, `testthat`:
  - de pre-flight zelftest (a)–(f) als unit-tests (schaalconventies
    `lavScores`, `information.observed`, α, T-tensor, V_inf);
  - regressietest HOIJ-2 vs. exacte bootstrap met *gedeelde* gewichten op een
    klein model (verschil = approximatiefout, moet klein en stabiel zijn);
  - equivalentietest IJ1 vs. gesloten-vorm invloedsfunctie.
- De **simulatiestudie zelf** (coverage, staartbalans, breedte, tijd vs.
  wald_inf/wald_hw/mc_hw/ij1/boot/BCa) is de wetenschappelijke validatie en
  bestaat al; het pakket moet dezelfde getallen reproduceren.

---

## 7. Stappenplan (prioriteit)

1. **Extraheer de §8-kern uit `HOIJ_simulatiestudie` naar een generieke
   functie** `hoijLavaan()` (§3): setup (S, H⁻¹, J, T) + gevectoriseerde
   replicatielus + SE/percentiel-CI + diagnostiek. Dit is vooral refactoren
   van bestaande, geteste code.
2. **Functionalen-interface**: `:=`-achtige expressies → φ(θ) op de vrije
   parameters (vervangt de hard-coded `functionals_scalar/vec`).
3. **Vervang `lavaan:::`-aanroepen** door `lav_export_estimation()` /
   geëxporteerde API; α-kalibratie wordt diagnostiek (§2).
4. **Scope-checks + meanstructure** (§4).
5. **Analytische casewise `Jᵢ`** via de Δ-Jacobiaan als performance-upgrade (§5).
6. **Pakket + tests + reproductie van de simulatieresultaten** (§6).
7. Fase 2: multi-group, missing (FIML), restricties.

---

### Referenties
- Giordano, Stephenson, Liu, Jordan & Broderick (2019), *A Swiss Army
  Infinitesimal Jackknife* (AISTATS) — IJ1; hogere-orde uitbreiding idem
  (Giordano et al., higher-order IJ).
- lavaan-internals: `lav_model_gradient.R`
  (<https://rdrr.io/cran/lavaan/src/R/lav_model_gradient.R>),
  `lav_model_vcov.R`
  (<https://github.com/yrosseel/lavaan/blob/master/R/lav_model_vcov.R>),
  `lav_export_estimation` (lavaan ≥ 0.6-17).
