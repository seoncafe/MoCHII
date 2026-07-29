# MoCHII plane-parallel 및 He I 10830 계획 검토 보고서

검토일: 2026-07-24  
검토 대상:

- `../MoCHII/docs/PLANE_PARALLEL_PLAN.md`
- `../MoCHII/docs/HEI_10833_PLAN.md`
- MoCHII commit `ff8a4e8`의 관련 구현

## 1. 요약

Plane-parallel의 핵심 수송은 대체로 올바르게 구현되어 있다.

- AMR과 Cartesian DDA 모두 x/y periodic wrap을 지원한다.
- AMR의 exact-face 재진입은 부동소수점 오차에 의한 photon loss를 방지한다.
- Beam의 `F_z = mu F_n` 및 isotropic field의 `F_z = pi I` 정규화가 올바르다.
- 직접광과 dust-scattered escape가 boundary `I(mu)` tally에 포함된다.
- 1-rank smoke test에서 `theta = 60 deg` beam의 전송률은 `0.78746`으로 기존 analytic 결과와 일치했다.

그러나 grazing ray 무한 루프, slab과 internal source의 합성, FUV 정규화의 의미, 양면 조명의 출력 표기 등은 추가 수정이 필요하다. 보수적 산란에 대한 Chandrasekhar benchmark와 자동 regression test도 아직 없다.

He I 10830 계획은 현재 형태로 구현하면 안 된다. 현재 사용 중인 Porter case-B He I 표는 단순 재결합 표가 아니라 metastable `2^3S`로부터의 collisional excitation을 이미 포함한 collisional-radiative emissivity이다. 따라서 현재 Porter 10830 출력에 별도의 `hec` 충돌 성분을 더하면 이중 계산이 된다.

He I 업데이트의 목표는 “누락된 충돌 방출 추가”가 아니라, Porter 표의 적용 범위를 벗어나는 환경에서 metastable `2^3S` population, column density, photoionization 및 10830 absorption을 명시적으로 계산하는 것으로 다시 정의하는 것이 적절하다.

## 2. Plane-parallel 구현 검토

### 2.1 올바르게 구현된 부분

#### Ionizing-band periodic transport

다음 경로에 x/y periodic wrapping이 적용되어 있다.

- `raytrace_ion_to_edge_amr`
- `raytrace_ion_to_tau_amr`
- `raytrace_ion_tau_only_amr`
- 대응하는 Cartesian DDA 경로

AMR 경로는 `slab_wrap_xy`에서 나간 좌표를 단순히 `x +/- xrange`로 이동하지 않고 반대쪽 경계의 정확한 `xmin`, `xmax`, `ymin`, `ymax` 값으로 재설정한다. 이는 경계 계산의 반올림 오차로 좌표가 box 밖에 놓여 `amr_find_leaf`가 packet을 버리는 문제를 방지한다.

관련 코드:

- `../MoCHII/src/raytrace_amr.f90:302`
- `../MoCHII/src/raytrace_amr.f90:313`
- `../MoCHII/src/raytrace_amr.f90:383`
- `../MoCHII/src/raytrace_amr.f90:491`
- `../MoCHII/src/raytrace_amr.f90:547`

#### Slab source normalization

`slab_setup`은 다음과 같이 수평면을 통과하는 flux를 계산한다.

```text
beam:      F_z = mu F_n,  mu = cos(theta)
isotropic: F_z = pi I
```

한 periodic tile에 들어오는 Monte Carlo luminosity는

```text
L_slab = F_z A_zface
```

로 설정된다. 이 정의는 물리적으로 올바르며, tile 면적이 packet weight와 packet 수의 면적 점유율 사이에서 상쇄되도록 한다.

관련 코드:

- `../MoCHII/src/ion_band_mod.f90:1640`
- `../MoCHII/src/ion_band_mod.f90:1664`
- `../MoCHII/src/ion_band_mod.f90:1673`

#### Boundary intensity tally

z boundary로 탈출하는 packet은 `mu = abs(kz)`에 따라 binning되며, 다음 식으로 azimuth-averaged intensity로 변환된다.

```text
I(mu) = L_escape_bin / (A_zface 2 pi mu Delta_mu)
```

따라서 수치 적분하면

```text
F = 2 pi integral I(mu) mu dmu
```

를 복원한다. 직접광과 scattered flight의 탈출이 모두 tally에 포함되어 있다.

관련 코드:

- `../MoCHII/src/jtally_mod.f90:68`
- `../MoCHII/src/jtally_mod.f90:95`
- `../MoCHII/src/raytrace_amr.f90:372`
- `../MoCHII/src/raytrace_amr.f90:437`
- `../MoCHII/src/raytrace_amr.f90:539`
- `../MoCHII/src/raytrace_amr.f90:606`

### 2.2 반드시 수정할 사항

#### P1. Grazing ray의 무한 루프

현재 lit beam face의 incidence angle 범위가 검증되지 않는다. `theta = 90 deg`이면

```text
mu = cos(theta) = 0
kz = 0
```

이 된다. 이 packet은 x/y 경계를 계속 wrap하지만 z 경계에는 도달하지 않는다. 경로 전체에서 opacity가 0이면 누적 optical depth도 증가하지 않으므로 AMR 및 DDA edge walk가 종료되지 않는다.

Isotropic sampling에서도 난수가 정확히 0인 경우 `mu = sqrt(u) = 0`이 될 수 있으므로 source setup 검증만으로 완전히 해결되지 않는다.

필요한 수정:

- 켜진 beam face마다 `0 <= theta < 90` 검사
- transport 내부에서 `abs(kz)`가 너무 작은 packet에 대한 명시적 처리
- 최대 periodic wrap 수 또는 최대 transport step 수 설정
- opacity가 0인 periodic grazing path를 종료시키는 guard

관련 코드:

- `../MoCHII/src/setup.f90:193`
- `../MoCHII/src/ion_band_mod.f90:1668`
- `../MoCHII/src/ion_band_mod.f90:1724`
- `../MoCHII/src/ion_band_mod.f90:1764`
- `../MoCHII/src/raytrace_amr.f90:413`
- `../MoCHII/src/raytrace_amr.f90:510`

#### P2. Slab과 internal source의 합성이 구현되지 않음

`source_geometry = 'slab'`을 선택하면 setup에서 `par%nsource = 0`으로 강제된다. 따라서 계획의 목표였던 externally illuminated slab 안의 embedded source를 한 run에서 사용할 수 없다.

또한 `source_geometry = 'slab'`과 기존 `ext_intensity`를 동시에 지정해도 ion source setup은 slab을 단일 component로 취급하므로 기존 external field가 사실상 무시된다. 현재 setup은 이를 오류나 경고로 알리지 않는다.

필요한 수정:

- `source_geometry`의 상호 배타적 분기 대신 `add_slab`과 같은 독립 source component로 구현하거나,
- slab face들을 기존 multi-component CDF에 직접 추가
- 지원 전까지는 slab과 기존 point/external source가 동시에 지정될 때 명시적 오류 또는 경고 출력

관련 코드:

- `../MoCHII/src/setup.f90:193`
- `../MoCHII/src/setup.f90:213`
- `../MoCHII/src/ion_band_mod.f90:167`
- `../MoCHII/src/ion_band_mod.f90:269`

#### P3. Internal-source periodic run에 boundary observable이 없음

`ion_peel`은 모든 `xy_periodic` run에서 금지된다. 그러나 boundary `I(mu)` tally는 `source_geometry = 'slab'`일 때만 활성화된다.

따라서 internal point source + `xy_periodic` 조합은:

- point-observer peel-off를 사용할 수 없고,
- boundary `I(mu)`도 출력되지 않는다.

이는 periodic slab의 observable을 boundary angular tally로 대체한다는 계획과 맞지 않는다.

권장 수정:

- `source_geometry`와 무관하게 `xy_periodic`이면 boundary tally를 활성화
- 입력 luminosity와 출력 header는 source configuration에 맞춰 일반화

관련 코드:

- `../MoCHII/src/setup.f90:379`
- `../MoCHII/src/main.f90:87`

#### P4. `add_fuv` 사용 시 `F_z`의 의미가 문서와 다름

Slab source에서는 `F_z A_zface`를 ionizing segment의 normalization으로 사용한다. `add_fuv`가 켜져 있으면 같은 spectrum의 FUV luminosity가 추가되어

```text
ion_Ltot = L_ion + L_FUV
```

가 된다. 따라서 실제 transported total flux는

```text
F_total = F_z (1 + L_FUV/L_ion)
```

이다.

그러나 계획 및 사용자 문서는 `F_z`를 전체 physics를 정하는 유일한 입사 flux로 설명한다. 특히 He metastable photoionization을 위해 4.8--13.6 eV FUV를 추가할 경우 이 차이가 중요해진다.

다음 중 하나로 의미를 확정해야 한다.

1. `slab_*_source`를 ionizing-band source strength로 정의하고 `F_ion`, `F_FUV`, `F_total`을 별도로 출력한다.
2. 전체 transported band가 지정된 `F_z`가 되도록 ionizing+FUV 전체 spectrum을 정규화한다.

관련 코드:

- `../MoCHII/src/ion_band_mod.f90:177`
- `../MoCHII/src/ion_band_mod.f90:182`
- `../MoCHII/src/ion_band_mod.f90:200`
- `../MoCHII/src/ion_band_mod.f90:540`

#### P5. Bottom/both illumination에서 reflection/transmission 표기가 잘못됨

출력 header는 항상

```text
L_top = reflected
L_bot = transmitted
```

로 기록한다.

- top-only illumination에서는 맞다.
- bottom-only에서는 reflection과 transmission의 의미가 반대다.
- both illumination에서는 top escape가 top 입사광의 reflection과 bottom 입사광의 transmission을 동시에 포함하므로 두 물리량을 분리할 수 없다.

권장 수정:

- 항상 유효한 기본 이름을 `L_escape_top`, `L_escape_bottom`으로 변경
- single-face illumination에서만 reflection/transmission alias 및 fraction 출력
- both-face run에서 양쪽 기원을 분리하려면 packet에 incident-face/component tag를 추가

관련 코드:

- `../MoCHII/src/jtally_mod.f90:104`
- `../MoCHII/src/jtally_mod.f90:109`

### 2.3 추가 권장 수정

#### 입력값 검증 강화

다음 setup 검사가 필요하다.

- `slab_nmu >= 1`
- lit face의 `slab_*_source > 0`
- lit beam face의 `0 <= theta < 90`
- `nx >= 1`, `ny >= 1`, `nz >= 1`
- `xmax`, `ymax`, `zmax > 0`

현재는 두 face의 합계 luminosity가 양수인지만 확인하므로 한 face의 음수 source가 다른 face의 양수 값에 가려질 수 있다. 이 경우 `slab_ptop`이 `[0,1]` 밖으로 나갈 수 있다.

#### 자동 regression test 추가

`tests/slab/`에는 입력 파일은 있지만 결과를 자동 판정하는 checker가 없다. 다음 gate를 스크립트로 고정할 필요가 있다.

- normal/oblique beam의 `exp(-tau/mu)`
- isotropic illumination의 `E_3(tau)`
- DDA와 shared/AMR 일치
- `nx = ny = 1`
- tile-area 및 lateral-resolution invariance
- bottom-only 및 both-face
- random/Sobol 일치
- MPI rank-count invariance
- conservative scattering의 Chandrasekhar benchmark
- grazing-ray fail-fast
- invalid `slab_nmu`, source 및 angle 검증

특히 `PLANE_PARALLEL_PLAN.md`에서도 conservative-scattering gate G-c는 아직 deferred 상태다.

#### 문서와 로그 정리

- 계획 문서 첫머리의 `Status: plan only, not started`는 실제 상태와 다르다.
- 계획의 `ext_geometry = 'slab'` 결정과 실제 `source_geometry = 'slab'` 구현을 통일해야 한다.
- slab run의 spectrum log가 다시 `external rec (isotropic)`이라고 표시되는 로깅 오류가 있다.
- `raytrace_ion_to_edge_amr`의 “No periodic boundaries” 주석은 현재 구현과 맞지 않는다.

## 3. He I 10830 계획 검토

### 3.1 핵심 전제 오류: Porter 10830은 이미 collisional excitation을 포함

`HEI_10833_PLAN.md`는 현재 MoCHII의 Porter 10830이 재결합 성분만 포함하고 collisional channel이 빠져 있다고 가정한다. 이 전제는 잘못되었다.

Porter et al. He I case-B emissivity는 collisional-radiative calculation이다. 특히 metastable `2^3S`에서 일어나는 collisional excitation이 10830 emissivity에 포함되어 있고, 그 결과 emissivity coefficient가 electron density에 강하게 의존한다.

저장소의 `hei_porter_caseB.txt`에서도 이를 직접 확인할 수 있다. `T = 10000 K`일 때 `4 pi j/(n_e n_He+)`는:

```text
n_e = 10 cm^-3:    2.8059e-25 erg cm^3 s^-1
n_e = 1e4 cm^-3:   1.8540e-24 erg cm^3 s^-1
```

로 약 6.6배 증가한다. 이 밀도 의존성은 표 안에 collisional excitation 및 metastable population 효과가 이미 들어 있음을 보여준다.

따라서 다음 방식은 이중 계산이다.

```text
current Porter HeI10830 + new hec collisional component
```

관련 파일:

- `../MoCHII/docs/HEI_10833_PLAN.md:4`
- `../MoCHII/data/atomic/hei_porter_caseB.txt:161`
- `../MoCHII/tools/fitting/make_hei_lines.py:1`
- `../MoCHII/src/sh95_mod.f90:1`
- `../MoCHII/src/lines_mod.f90:170`

참고 문헌:

- Porter et al. 2012, *Improved He I Emissivities in the Case B Approximation*: <https://arxiv.org/abs/1206.4115>
- Del Zanna & Storey 2022, *Helium line emissivities for nebular astrophysics*: <https://academic.oup.com/mnras/article/513/1/1198/6564189>

### 3.2 OH18 population 식의 수정

OH18은 다음 상태를 사용한다.

- neutral singlet ground state `f1`
- metastable triplet state `f3`
- He ion

계획에 적힌 “`1^1S`, `2^3S`, `2^1S` plus ion”의 독립 3-level neutral model은 OH18 설명과 다르다.

OH18의 triplet equation을 local steady-state로 단순화하면 다음과 같다.

```text
source:
    alpha_3 n_e n(He+)
  + q_13a n_e n(He 1^1S)

sinks:
  n_3 [
      A_31
    + Phi_3
    + n_e (q_31a + q_31b)
    + n_HI Q_31
  ]
```

계획에서 수정할 사항:

- OH18 식에는 `q_13b`가 없다.
- `q_13a`는 singlet ground에서 triplet으로 들어가는 source다.
- `q_31a + q_31b`가 electron-impact triplet sink다.
- H collision 항은 총 `n_H`가 아니라 neutral `n_HI`를 사용해야 한다.
- `Q_31`은 단순 de-excitation이 아니라 associative ionization과 Penning ionization의 합이다.

참고 문헌:

- Oklopčić & Hirata 2018, *A New Window into Escaping Exoplanet Atmospheres*: <https://arxiv.org/abs/1711.05269>

### 3.3 Local equilibrium과 exoplanet advection을 구분해야 함

OH18의 실제 식은 local algebraic equilibrium이 아니라 advection equation이다.

```text
v df_3/dr = source - sink
```

현재 계획처럼 cell마다 `source = sink`를 적용하면 static/local-equilibrium approximation이 된다. 이는 dense H II region이나 반응 시간척도가 flow 시간척도보다 충분히 짧은 조건에서는 사용할 수 있지만, OH18 hot-Jupiter benchmark를 그대로 재현한다고 볼 수 없다.

현재 MoCHII의 generic AMR reader는 입력 velocity column을 무시한다. 따라서 계획의 exoplanet G3를 수행하려면 다음 중 하나가 필요하다.

- velocity/streamline 입력 및 advection ODE 구현
- 또는 OH18 profile의 density, temperature, velocity를 읽어 별도 1D post-processing solver로 적분

Advection을 구현하지 않는다면 G3를 local-equilibrium approximation test로 이름과 범위를 바꿔야 한다.

관련 코드:

- `../MoCHII/src/read_generic_amr.f90:18`

### 3.4 OH18은 주로 absorption 모델

OH18의 주 목적은 metastable helium의 10830 transit absorption 계산이다. `2^3S <-> 2^3P` 전이는 excitation 후 거의 다시 `2^3S`로 돌아오므로 metastable reservoir를 유의미하게 소모하지 않아 population network에서 제외된다.

따라서 OH18 population model을 도입할 때 우선적인 observable은:

- `n_2s3`
- metastable column density
- 10830 absorption coefficient
- optical depth와 transit absorption

이다. 단순한 optically thin collisional emission row 추가는 exoplanet 목적과 직접 일치하지 않는다.

### 3.5 FUV photoionization은 exoplanet 적용에서 필수

`2^3S` photoionization threshold는 약 4.77 eV다. Dense H II region을 위한 첫 근사에서는 `Phi_3 = 0`을 사용할 수 있지만, exoplanet metastable population은 stellar near-UV/FUV field에 매우 민감하다.

따라서:

- dense nebular mode에서는 `Phi_3 = 0` 근사를 허용할 수 있다.
- exoplanet/PDR mode에서는 `Phi_3`를 optional 후속 단계로 두면 안 된다.
- transported `J_nu`를 이용해 다음 적분을 계산해야 한다.

```text
Phi_3 = 4 pi integral J_nu sigma_3(nu)/(h nu) dnu
```

이때 `efuv_min`뿐 아니라 bin edge가 4.77 eV threshold를 실제로 포함하는지 확인해야 한다.

### 3.6 CHIANTI pipeline만으로는 충분하지 않음

CHIANTI `he_1`의 `elvlc`, `wgfa`, `scups`는 bound-level energy, radiative transition 및 electron collision 자료를 제공할 수 있다. 그러나 explicit metastable model에는 추가로 다음 자료가 필요하다.

- triplet recombination coefficient `alpha_3(T)`
- metastable photoionization cross section
- neutral-H associative/Penning ionization `Q_31(T)`
- 필요할 경우 collisionless recombination contribution to 10830

따라서 일반 `nlevel_he_1.txt` 생성만으로 population과 total 10830 emissivity를 일관되게 계산할 수 없다.

## 4. 권장 He I 설계

### 4.1 모델을 두 가지로 분리

#### `porter_caseB`

현재 Porter table을 그대로 사용한다.

- H II region용 case-B collisional-radiative total emissivity
- density 및 temperature dependence 포함
- 기본값으로 유지

문서와 코드 주석의 “recombination-only” 표현은 “case-B collisional-radiative emissivity”로 수정한다.

#### `explicit_meta`

환경별 metastable population을 별도로 계산한다.

출력:

- `n_2s3`
- `f_2s3`
- `Phi_3`
- metastable column
- 10830 absorption opacity 또는 optical depth

Local equilibrium과 advection mode는 별도로 구분한다.

### 4.2 Explicit emission은 Porter row를 대체해야 함

Explicit metastable emission을 제공하려면 다음과 같은 분리가 필요하다.

```text
j_total(10830)
  = j_recomb_collisionless(10830)
  + n_3 n_e q(2^3S -> 2^3P) h nu B
```

여기서 `j_recomb_collisionless`는 현재 Porter total table과 다른 자료여야 한다. 현재 Porter coefficient를 첫 항으로 사용하면 이미 포함된 collisional excitation을 다시 더하게 된다.

권장 동작:

- `porter_caseB` mode에서는 현재 `HeI10830` row만 출력
- `explicit_meta` mode에서는 Porter `HeI10830` row를 대체
- 두 결과를 동시에 합산하지 않음
- 필요하면 diagnostic 목적으로 `recomb`과 `coll`을 별도 column으로 출력하되 total 정의를 명확히 함

### 4.3 수정된 validation gates

#### H1. Current Porter baseline

- 현재 MoCHII interpolation이 원본 Cloudy/Porter grid를 정확히 재현하는지 검사
- 여러 `T_e`, `n_e`에서 10830 및 다른 He I line 비교

#### H2. Explicit local-equilibrium population

- OH18 rate coefficients를 사용한 one-zone algebraic solution과 비교
- `Phi_3`, electron collision, neutral-H destruction을 각각 켜고 끄는 limit test

#### H3. Porter-equivalent limit

- Porter가 가정한 case-B 조건에 가까운 환경에서 explicit total emissivity가 Porter table과 일치하는지 검사
- 기존 계획처럼 “Porter + explicit collisional”을 Cloudy와 비교하지 않음

#### H4. Exoplanet advection

- velocity/advection solver를 구현한 뒤 OH18 hot-Jupiter profile 및 metastable column과 비교
- local-equilibrium mode만 구현한 경우에는 G3 통과를 주장하지 않음

#### H5. Absorption and line transfer

- optically thin column-to-optical-depth analytic test
- 10830 fine-structure triplet oscillator strengths 포함
- line-center optical depth가 큰 경우 단순 emissivity/column 결과에 경고

## 5. 권장 작업 순서

### Plane-parallel

1. Slab parameter sanity checks와 grazing guard 추가
2. `L_escape_top/bottom` 출력 의미 수정
3. FUV 포함 시 flux normalization convention 확정
4. 모든 `xy_periodic` source에 boundary tally 제공
5. Slab을 multi-source component로 일반화
6. 자동 regression checker 및 conservative-scattering benchmark 추가

### He I 10830

1. 현재 Porter 표가 collisional-radiative total임을 문서와 코드에 명시
2. 기존 `HEI_10833_PLAN.md`의 “missing collisional component” 전제 제거
3. `hei_meta_mod`를 emission addition이 아닌 metastable diagnostic으로 설계
4. OH18 source/sink 식과 neutral-H 항 수정
5. `Phi_3`를 포함한 local-equilibrium population 구현
6. `n_2s3`, column, opacity 출력
7. Explicit emission이 필요하면 collisionless recombination 자료를 확보하여 Porter 10830 row를 대체
8. Exoplanet benchmark가 필요하면 velocity/advection solver 추가

## 6. 최종 판단

Plane-parallel 구현은 top-only 또는 bottom-only의 단일 slab illumination, 유한한 `abs(kz)`, 현재 band normalization이라는 제한 아래에서는 핵심 알고리즘이 올바르다. 다만 현재 계획 전체의 목표인 arbitrary source composition과 모든 periodic slab의 boundary observable까지는 완성되지 않았다.

He I 10830 계획은 핵심 전제가 잘못되어 전면 수정이 필요하다. 현재 Porter 10830에 별도 collisional row를 더해서는 안 된다. 새 기능의 핵심 가치는 Porter 표에 이미 포함된 nebular collisional emission을 반복 계산하는 것이 아니라, FUV radiation 및 neutral-H destruction이 중요한 환경에서 metastable helium population과 10830 absorption을 물리적으로 계산하는 데 있다.
