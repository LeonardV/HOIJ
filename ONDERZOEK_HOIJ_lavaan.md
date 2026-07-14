# Onderzoek: HOIJ implementeren in lavaan

Dit document beschrijft wat er nodig is om de HOIJ-procedure als volwaardige
implementatie boven op `lavaan` te bouwen. Het vertrekt vanuit de bestaande
prototype-functies in `kernfuncties` en brengt in kaart wat af is, wat ontbreekt,
en welke keuzes/risico's er liggen.

> Terminologie: het bestand `kernfuncties` levert de *rekenkern* (casewise
> log-likelihood, casewise informatiebijdragen `J`, de derde-orde tensor `T`,
> en een kalibratiestap). Dit zijn precies de bouwstenen van een **hogere-orde
> correctie**: de gebruikelijke informatie­matrices (verwacht `I` / geobserveerd
> `J`) aangevuld met de derde-afgeleide tensor `T`. De uiteindelijke HOIJ-grootheid
> (de gecorrigeerde variantie / toetsgrootheid / betrouwbaarheidsinterval) zelf
> staat nog niet in het bestand — dat is het grootste inhoudelijke gat (zie §2).

---

## 1. Wat de huidige rekenkern doet (analyse van `kernfuncties`)

| Functie | Wat het berekent | Schaal / conventie |
|---|---|---|
| `compute_loglik_casewise(fit, theta)` | Per observatie de multivariaat-normale log-likelihood ℓᵢ(θ). Zet `theta` (lavaan's vrije-parametervector `x`) om via `lav_model_x2glist` → model-implied Σ(θ), μ(θ), en evalueert de MVN-dichtheid per rij. Geeft vector van lengte N. | log-likelihood |
| `compute_all_J(fit, theta0, delta)` | N×D×D array met `J[i,,] = −∂²ℓᵢ/∂θ∂θᵀ`: de casewise **geobserveerde informatie**. Som over i = totale geobserveerde informatiematrix. Centrale tweede finite differences. | log-likelihood |
| `make_grad_F(fit)` | Closure `grad_F(theta) = ∂F/∂θ`, met lavaan's **analytische** gradiënt van de discrepantie-/fitfunctie F. | discrepantie F |
| `calibrate_alpha(grad_F, theta0, H_observed, h)` | Differentieert `grad_F` numeriek → H_grad = ∂²F/∂θ∂θᵀ, en zoekt de scalar α met `H_observed ≈ α·H_grad`. `spread` = hoe constant de ratio is (consistentie-check). | brug tussen beide schalen |
| `compute_T_tensor_grad(grad_F, theta, alpha, h)` | D×D×D tensor `T = α·∂³F/∂θ∂θ∂θ`, via tweemaal numeriek differentiëren van `grad_F`, daarna volledig gesymmetriseerd over alle 6 permutaties. | log-likelihood (via α) |

**Kernobservatie over de schalen.** `J` staat op de log-likelihood-schaal,
`grad_F` op de discrepantie-schaal F (voor ML geldt ruwweg F ≈ −(1/N)·Σℓᵢ + const).
`calibrate_alpha` overbrugt die twee schalen empirisch (verwacht α ≈ N), en `T`
wordt met diezelfde α teruggeschaald naar de log-likelihood-schaal. Deze
schaalbrug is nu *empirisch gekalibreerd*; voor een productie-implementatie zou
alles op één consistente (log-likelihood) schaal analytisch berekend moeten
worden, zodat α niet meer nodig is behalve als diagnostiek.

---

## 2. Grootste inhoudelijke gat: de eigenlijke HOIJ-grootheid ontbreekt

De rekenkern levert `J` (geobserveerde informatie per case) en `T` (derde-orde
tensor), maar **nergens worden die samengevoegd tot het eindresultaat**. Nodig:

- Een expliciete definitie + implementatie van de doelgrootheid: de hogere-orde
  gecorrigeerde (co)variantie van θ̂, of de gecorrigeerde toetsgrootheid /
  het gecorrigeerde betrouwbaarheidsinterval (bv. Edgeworth-/Cornish-Fisher- of
  Bartlett-achtige correctie die `T` gebruikt om scheefheid/bias te corrigeren).
- De **verwachte informatie `I`** ontbreekt (er is alleen geobserveerde `J`).
  Als HOIJ `I` tegen `J` afzet, is die apart nodig — via
  `lavInspect(fit, "information.expected")` of `lavaan:::lav_model_information`.
- Eventueel de **casewise scores** (outer product = "meat" van een sandwich):
  officieel beschikbaar via `lavaan::lavScores(fit)` / `estfun`, in plaats van
  numeriek afgeleid.

Dit is de eerste prioriteit: zonder de assemblagestap is er nog geen HOIJ-output.

---

## 3. Afhankelijkheid van niet-geëxporteerde lavaan-internals (risico)

De kern leunt op `lavaan:::`-functies: `lav_model_x2glist`, `lav_model_implied`,
`lav_model_gradient`. Dat zijn interne functies die tussen lavaan-versies kunnen
veranderen (signatuur, gedrag, slot-namen). Aanbevelingen:

- Gebruik waar mogelijk de **officieel geëxporteerde** hooks. `lavaan` biedt
  `lav_export_estimation(fit)` dat objective- en gradiëntfuncties naar buiten
  geeft; gebruik dat i.p.v. `make_grad_F` met `:::`.
- Casewise scores en informatie via geëxporteerde API: `lavScores()`,
  `lavInspect(fit, "information")`, `lavInspect(fit, "information.observed")`,
  `lavInspect(fit, "information.expected")`, `vcov()`.
- Waar `:::` onvermijdelijk blijft: pin een minimale lavaan-versie, voeg een
  versie-check toe, en dek af met tests die breken als de internals wijzigen.

---

## 4. Dekkingsgaten t.o.v. lavaan's generaliteit (grote implementatie-oppervlak)

De prototype-kern maakt sterke, impliciete aannames. Voor elke lavaan-situatie
die je wil ondersteunen moet de kern uitgebreid worden:

- **Eén groep.** Overal wordt `[[1]]` gebruikt (`X[[1]]`, `cov[[1]]`, `mean[[1]]`).
  Multi-group vereist een lus over groepen en blok-combinatie van J/T.
- **Volledige data + ML.** Geen FIML/missing (`missing="ml"`), geen
  categorisch/ordinaal (WLSMV/DWLS), geen niet-normaal-robuust (MLR/MLM), geen
  sampling weights, geen clustered/multilevel. Elk daarvan verandert de casewise
  likelihood én de afgeleidenstructuur. Beslis de scope en **geef een nette fout**
  bij niet-ondersteunde estimators i.p.v. stilzwijgend fout te rekenen.
- **Meanstructure.** Nu ad-hoc: bij ontbrekende μ valt de code terug op
  `colMeans(X)`. Moet netjes `meanstructure`, `fixed.x`, `conditional.x` en
  intercepts respecteren.
- **Constraints en gedefinieerde (`:=`) parameters.** `theta` is hier de vrije
  `x`-vector. Bij (on)gelijkheidsrestricties moet de derde-orde tensor op de
  *gerestricteerde* variëteit worden uitgedrukt (constraint-Jacobiaan
  `lavmodel@con.jac`), of moet je herparametriseren. Dit is precies het raakvlak
  met restriktor/inequality-constrained inferentie en verdient expliciete aandacht.

---

## 5. Numerieke robuustheid

- **Finite differences zijn fragiel.** De tweede- (`delta²`) en derde-orde
  (`h²·…`) schema's versterken ruis; stapgroottes staan hard-coded. Overweeg
  Richardson-extrapolatie (`numDeriv::genD`/`hessian`) of, beter,
  **analytische afgeleiden** waar lavaan die al levert (Δ, de delta-/Jacobiaan
  van de implied moments). Analytische casewise scores en geobserveerde
  informatie zijn nauwkeuriger én sneller.
- **Kosten.** De `J`-lus herberekent de volledige casewise loglik O(D²) keer;
  `T` roept `grad_F` O(D²) keer aan. Voor realistische D is dat traag. Analytische
  afgeleiden of vectorisatie zijn nodig voor bruikbare snelheid.
- **Randgevallen** die afgevangen moeten worden: niet-positief-definiete Σ, bijna
  singuliere Hessiaan (nu al `MASS::ginv`-fallback), parameters op de rand,
  slechte kalibratie (`spread` groot → waarschuwen/afbreken).

---

## 6. Software-engineering om er een lavaan-waardig gereedschap van te maken

- **R-pakketstructuur** (standalone of restriktor-adjacent): `DESCRIPTION`,
  `NAMESPACE`, `Imports: lavaan`, documentatie (roxygen), `testthat`-tests.
- **Publieke API**: één functie die een *gefitte* lavaan-object neemt, de
  aannames vooraf valideert (estimator, data, groepen, constraints), en de
  HOIJ-uitvoer teruggeeft in een net object met `print`/`summary`.
- **Numerieke safeguards** en informatieve warnings/errors i.p.v. stille NA's.

---

## 7. Validatie & verificatie

- Behoud en veralgemeen de `calibrate_alpha`-`spread` als interne
  consistentie-check (numerieke vs. analytische Hessiaan; verwacht α ≈ N).
- Cross-check de numerieke derde afgeleide tegen een symbolische/analytische
  referentie op een klein model; controleer symmetrie/invariantie van `T`.
- **Simulatiestudie**: toon dat de hogere-orde correctie daadwerkelijk de
  dekking (coverage) / type-I-fout verbetert t.o.v. de standaard Wald/LR, en
  benchmark tegen bootstrap.

---

## 8. Concreet stappenplan (prioriteit)

1. **Definieer en implementeer de assemblagestap** (§2): de eigenlijke HOIJ-formule
   die `J`, `I` en `T` combineert tot de gecorrigeerde variantie/toets. Zonder dit
   is er geen output.
2. **Stap over op geëxporteerde lavaan-API** (`lav_export_estimation`, `lavScores`,
   `lavInspect(..., "information*")`) om de `:::`-afhankelijkheid en de empirische
   α-kalibratie te vervangen door één consistente, analytische schaal (§3, §1).
3. **Vervang finite differences door analytische/Richardson-afgeleiden** voor
   nauwkeurigheid en snelheid (§5).
4. **Bepaal de scope en bewaak die**: begin met één groep + complete data + ML,
   met harde checks die niet-ondersteunde estimators netjes afwijzen (§4).
5. **Constraints/`:=`-parameters** correct afhandelen op de gerestricteerde
   variëteit (§4) — sluit aan bij het restriktor-werk.
6. **Verpak als R-pakket met tests + simulatie-validatie** (§6, §7).

---

### Bronnen (lavaan-internals & informatiematrix-machinerie)
- lavaan model-functies: <https://rdrr.io/cran/lavaan/man/lav_model.html>
- `lav_model_gradient.R`: <https://rdrr.io/cran/lavaan/src/R/lav_model_gradient.R>
- `lav_model_vcov.R` (informatie/vcov/sandwich): <https://github.com/yrosseel/lavaan/blob/master/R/lav_model_vcov.R>
- lavaan-class (slots): <https://rdrr.io/cran/lavaan/man/lavaan-class.html>
