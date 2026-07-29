# MoCHII 코드 감사 및 부동소수점 분석

- 작성일: 2026-07-13
- 대상: `c54d270` 기준 소스와 본 감사에서 적용한 수정
- 컴파일러: Intel Fortran 2025.3.2 (`mpiifort`/`mpiifx`, IFX 계열)
- 주요 실행 검증: 8 MPI ranks, `tests/peel/peel_dbg.in`, rank 합계 2,000,000 photon

## 1. 결론 요약

`ion_peel_mod.f90`만 `-fp-model precise`가 필요했던 직접 원인은 peel 누적 자체가
아니라 다음 두 수치 결함의 조합이었다.

1. `raytrace_amr.f90`의 극축 산란 기저가 직교하지 않아 산란을 거친 방향벡터의
   norm이 1에서 벗어났다.
2. HG phase function 분모의 보호식 `max(base, tiny())**1.5`에서
   `tiny()**1.5`가 double precision 범위 아래로 underflow하여 0이 되었다.

과거 재현 로그에는 실제로 `cos(theta)=1.00745`, `peel=Infinity`가 남아 있다
([`tests/peel/peel_dbg.log`](../tests/peel/peel_dbg.log), PEEL 부분). 이는 단순한
마지막 자리 반올림이 아니라 단위 방향벡터 불변식이 깨졌다는 증거다.

두 경로를 수정한 뒤 `ion_peel_mod.o`의 `precise` 예외를 Makefile에서 제거했다.
수정된 코드의 `fast=2`와 `precise` 결과는 모든 image pixel이 finite였고, 네 채널
합계의 상대 차이가 `5.7e-16`에서 `1.0e-15`였다. 따라서 현재 검증 범위에서는
`ion_peel_mod.f90`도 다른 routine과 동일하게 `-fp-model fast=2`로 빌드할 수 있다.

단, 이는 `fast=2`가 일반적으로 과학 계산에 안전하다는 뜻은 아니다. Intel은
`fast=2`를 value safety가 없는 "very unsafe" 모드로 분류한다
([Intel Floating-Point Optimizations](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-1/floating-point-optimizations.html)).
전체 회귀 테스트를 컴파일러/CPU별로 유지해야 한다.

`~/CLOUDY/c23.01`의 atomic data와 추가 비교한 결과, photoionization cross
section과 CHIANTI line data는 즉시 교체할 이유가 없었다. 반면 고차
charge-transfer 계수에는 이온 stage를 한 칸 잘못 대응시킨 항목이 다수 있으며,
이는 atomic data 갱신보다 우선해서 수정해야 한다. Huang et al. (2023) 원문과
`EXHALE`의 63개 계수/반응 stage에서는 전사 오류를 찾지 못했고, 문제는 MoCHII가
EXHALE 계수에 MOCASSIN 고차 계수를 추가하는 과정에서 생겼다. 상세 내용은 9절에
정리했다.

수소 recombination line도 Cloudy와 계산 경로가 다르다. MoCHII는 모든 cell에서
Storey & Hummer (1995, SH95) case-B emissivity를 직접 보간하지만, Cloudy의 주
`H  1` Balmer line은 level-resolved H-like model atom의 population solve 결과다.
같은 SH95 reference인 Cloudy `Ca B`와 MoCHII 표는 0.2% 이내로 일치하므로 표 자체의
갱신보다 case-A/B 연동, 유효 범위 검사와 collisional H I 성분이 우선이다. 상세
내용은 10절에 정리했다.

## 2. 코드 구조와 분석 범위

주 실행 흐름은 다음과 같다.

1. 입력 및 MPI/shared-memory 초기화: `setup.f90`, `memory_mod_mpi.f90`
2. AMR 또는 uniform Cartesian grid 생성: `octree_mod.f90`, `grid_mod_amr.f90`
3. ionizing/FUV band, gas/dust opacity 구성: `ion_band_mod.f90`, `gas_opacity_mod.f90`
4. photon transport와 path-length tally: `raytrace_amr.f90`, `jtally_mod.f90`
5. H/He/metal ionization 및 thermal balance: `ion_balance_mod.f90`,
   `species_mod.f90`, `thermal_mod.f90`, `cooling_mod.f90`
6. line/continuum/dust/image 출력: `lines_mod.f90`, `nebcont_mod.f90`,
   `sedust_mod.f90`, `ion_peel_mod.f90`

수행한 점검은 다음과 같다.

- 약 2만 줄의 Fortran 소스에 대한 정적 검토와 수치 위험 패턴 검색
- Intel `-warn all -stand f18 -check all -fpe0 -O0` 전체 빌드
- 과거 peel 실패 로그와 git 변경 이력 검토
- `ion_peel_mod`/`raytrace_amr`의 fast/precise 분리 빌드와 MPI 실행 비교
- 기존 테스트 구조와 자동 gate 범위 검토
- Huang et al. (2023) 로컬 PDF Table 4와 EXHALE의 63개 charge-exchange 식,
  donor/acceptor stage 및 MoCHII/MOCASSIN 이식 경로 대조
- Cloudy H-like level population/emissivity 호출 경로와 MoCHII SH95 case-B
  보간의 직접 비교

전체 물리 benchmark를 모두 재실행하지는 않았다. 이 보고서의 물리 정확도 평가는
기존 test gate의 존재와 구현 검토에 기반하며, 신규 실행 비교는 peel debug case에
집중했다.

## 3. `fast=2` 오류의 상세 원인

### 3.1 비직교 산란 기저

기존 코드는 입사 방향 `k`에 수직인 벡터 `u`를 만들 때 `|kz| >= 0.99`이면
`u=(1,0,0)`을 사용했다. 일반적인 `k=(kx,ky,kz)`에서 `k dot u = kx`이므로 이
벡터는 수직이 아니다. 이 상태에서 `v=k cross u`와 HG 산란각을 결합하면 출력
방향벡터도 단위벡터가 아니며, 반복 산란 중 오차가 누적된다.

수정 내용은 [`src/raytrace_amr.f90`](../src/raytrace_amr.f90)의 HG 산란부에 있다.

- 입력 `k`를 먼저 정규화
- 극축에서 `u=(0,-kz,ky)`를 사용하여 `k dot u = 0` 보장
- 산란 후 출력 방향을 다시 정규화

이는 이미지 방어 코드가 아니라 transport 자체의 물리 정확성 수정이다. 기존
방향 norm 오류는 peel image뿐 아니라 이후 cell crossing 거리와 scattered
path-length tally에도 영향을 줄 수 있었다.

### 3.2 HG 분모 보호식의 underflow

HG phase function은

```text
P(cos theta) = (1-g^2) / [4 pi (1+g^2-2g cos theta)^(3/2)]
```

이다. 기존 `max(denominator, tiny())**1.5`는 `tiny()` 자체는 양수이지만 이를
`1.5` 제곱할 때 0으로 underflow하므로 division by zero를 막지 못한다.

수정된 [`src/ion_peel_mod.f90`](../src/ion_peel_mod.f90)은 다음을 적용한다.

- 입사 방향 norm을 반영하여 실제 cosine을 계산
- cosine을 `[-1,1]`로 clamp
- `g`를 `[-0.9999,0.9999]`로 clamp
- 분모의 물리적 최솟값 `(1-|g|)^2`를 사용
- `den**1.5` 대신 `den*sqrt(den)` 사용

마지막 변경은 underflow 원인을 제거하고 compiler가 분수 지수 연산을 변형할
여지도 줄인다.

### 3.3 fast-math와 NaN 검사

IFX에서 `fast`/`fast=2`는 NaN이 없다고 가정하는 최적화를 허용한다. Intel은
IFX에서 fast mode와 NaN 비교를 함께 사용할 때 `-assume nan_compare` 사용을
안내한다
([Intel Fortran Get Started Guide](https://www.intel.com/content/www/us/en/docs/fortran-compiler/get-started-guide/2025-0/overview.html)).
현재 peel의 `nonfinite()`는 실수 비교 대신 exponent bit를 검사하므로 이 변환의
영향을 받지 않는다.

반면 `utility.f90`의 `is_finite32/64()`는 `ieee_class`를 사용하며, 여러 입력의
기본값은 NaN sentinel이다. 모든 파일을 계속 `fast=2`로 컴파일하려면 다음 중
하나를 권장한다.

1. `-assume nan_compare`를 Intel release flag에 추가하고 회귀 테스트한다.
2. 더 확실한 방법으로 `utility.f90`의 finite 판정도 real32/real64 exponent-bit
   검사로 바꾼다.
3. 이 작은 utility만 `precise`로 두되, hot transport routine은 모두 fast로 둔다.

이번 peel 문제에는 2번 방식이 이미 국소적으로 적용되어 있으며, 별도 precise
object는 더 이상 필요하지 않다.

## 4. 실행 검증 결과

동일한 수정 소스를 다음 두 방식으로 컴파일해 비교했다.

```text
fast:    -ipo -O3 -no-prec-div -fp-model fast=2
precise: -ipo -O3 -fp-model precise
```

| image channel | fast 합계 | precise 합계 | 상대 차이 |
|---|---:|---:|---:|
| `direct_fuv` | `1.0027040687495079e+02` | `1.0027040687495074e+02` | `5.67e-16` |
| `direct_ion` | `4.3282976745793173e-48` | `4.3282976745793136e-48` | `8.42e-16` |
| `scatt_fuv` | `8.3178379379378086e+01` | `8.3178379379378015e+01` | `8.54e-16` |
| `scatt_ion` | `3.2987257096900135e-50` | `3.2987257096900102e-50` | `1.01e-15` |

- 두 실행 모두 exit code 0
- 두 실행 모두 non-finite pixel 0
- 전체 출력 로그의 peel 합계는 양쪽 모두 `direct=2.8357e29`, `scatt=2.3523e29`
- wall time은 fast 1.597분, precise 1.612분이었으나 단일 실행 차이이므로 성능
  결론으로 사용하면 안 된다.
- 모든 최종 수정 후 만든 fast 통합 build도 같은 8-rank case를 재실행해 exit code
  0, 네 channel 모두 finite, 위 fast channel 합계와 bit-for-bit 같은 결과를 확인했다.

## 5. 발견한 오류와 위험

### P0: 수정 완료

#### 산란 방향의 단위벡터 불변식 위반

위 3.1절의 문제다. scattered transport와 image 양쪽에 영향을 주는 실제 계산
오류였으며 수정했다.

#### HG 분모의 무효한 underflow guard

위 3.2절의 문제다. `precise`가 증상을 가렸지만 수식 자체를 안정화하는 것이 올바른
해결이다. 수정 후 Makefile의 object별 `precise` 옵션을 제거했다.

#### 비정방 3D density centering의 z축 인덱스 오타

[`src/read_mod.f90`](../src/read_mod.f90)의 `loc(3) < n3cen` branch가 z축 offset에
`n2cen`을 사용했다. `ny != nz`인 grid에서 잘못된 section을 복사하거나 bounds
오류를 낼 수 있다. `n3cen`으로 수정했다. 현재 test suite에는 이 branch를
검증하는 비정방 density cube test가 없다.

### P1: 수정 권장

#### 최대 iteration 도달 후에도 unconverged 결과를 정상 출력

[`src/main.f90`](../src/main.f90)의 iteration loop는 마지막 iteration까지
`converged=.false.`여도 경고 없이 imaging과 science output을 계속한다. 실제
`peel_dbg` 실행도 6회 뒤 `max delta Te/Te = 7.73`, volume 기준 `2.90e-2`인 상태로
출력했다. debug input의 짧은 iteration은 의도적일 수 있지만, production run에서
같은 상황을 구분할 metadata가 없다.

권장 조치:

- loop 종료 사유와 최종 convergence metric을 출력 파일 header에 기록
- 최대 iteration 도달 시 rank 0에서 명확한 warning 출력
- 선택적으로 `par%require_convergence`가 true이면 nonzero exit

#### 입력값 검증 부족과 알 수 없는 단위의 묵시적 `kpc` 처리

[`src/setup.f90`](../src/setup.f90)은 band/grid의 일부만 검증한다. 최소한
`no_photons>0`, `0<ion_relax<=1`, `0<te_min<te_max`, `nxim/nyim>0`,
`distance>0`, `case_ab in {A,B}`를 확인해야 한다. 특히 알 수 없는
`distance_unit`을 오류로 처리하지 않고 `kpc`로 바꾸는 현재 동작은 결과 단위를
조용히 `10^3` 이상 틀리게 만들 수 있다.

#### He equilibrium의 overflow 가능성

[`src/ion_balance_mod.f90`](../src/ion_balance_mod.f90)은
`xHeI=1/(1+r1+r1*r2)`, `xHeII=xHeI*r1`을 직접 계산한다. 매우 작은 `ne`와 큰
photo-rate 조합에서는 `r1*r2`가 overflow하고 이어지는 `0*Inf`가 NaN이 될 수
있다. scale을 나누어 쓰거나 log-ratio/정규화된 3-state weight 방식으로 바꾸는
것이 안전하다.

### P2: 견고성 개선

- `read_3D()`는 파일 open/read 실패 후에도 미할당 `arr`를 사용하려고 한다.
  오류 메시지만 출력하지 말고 즉시 `error stop` 또는 MPI abort해야 한다.
- 여러 MPI collective와 I/O call이 `ierr/status`를 확인하지 않는다. 실패한 output이
  정상 완료처럼 보일 수 있으므로 공통 error wrapper가 필요하다.
- explicit observer 좌표 mode에서 모든 observer에 첫 observer 기반의 공통
  `par%distance`를 사용한다. 서로 다른 반경의 observer 목록은 `cosb` 범위를
  벗어나 `sqrt(1-cosb^2)`를 invalid로 만들 수 있다. observer별 distance를 쓰거나
  입력 제약을 검사해야 한다.
- `utility.f90`의 NaN sentinel 판정은 fast-math/compiler 조합별 test가 필요하다.

## 6. 성능 개선 후보

### 1순위: cross section을 bin별로 사전 계산

[`src/gas_rates_mod.f90`](../src/gas_rates_mod.f90)은 모든 `(leaf,bin)`에서
`sigma_HI`, `sigma_HeI`, `sigma_HeII`를 다시 호출한다. cross section은 energy
bin에만 의존한다. 262,144 leaves와 48 bins이면 iteration마다 약 3,775만 회의
불필요한 함수 평가가 발생하며 내부에는 제곱근과 실수 지수가 있다.

`ion_setup()` 또는 `gas_opacity_setup()`에서 세 배열을 한 번 계산해 공유하면
정확도를 바꾸지 않고 큰 속도 향상을 기대할 수 있다. `1-Eth/E`, `1/(E*ev2erg)`도
같이 사전 계산할 수 있다.

### 2순위: equilibrium/thermal cell solve의 rank 중복 제거

`jt_ion`을 `MPI_ALLREDUCE`한 뒤 모든 MPI rank가 모든 leaf의 ion/thermal solve를
동일하게 반복하고, node의 `h_rank==0`만 shared state를 쓴다. 8 ranks에서 같은
262,144-cell bisection을 8번 수행하는 구조다.

단기 개선은 node마다 `h_rank==0` 하나만 solve하고 barrier하는 것이다. 장기적으로는
leaf를 rank에 분할하고 결과를 shared window에 직접 기록한 뒤 convergence scalar만
reduce하는 편이 확장성이 좋다.

### 3순위: tally 메모리와 collective의 node-aware화

검증 case에서 `jt_ion`은 rank당 0.094 GB였다. 32 ranks면 node당 약 3 GB이며,
매 iteration 전체 배열을 `MPI_ALLREDUCE`한다. shared-memory node reduction 후
동일 local-rank communicator로 inter-node reduction하는 기존 `memory_mod_mpi`
패턴을 tally에도 적용하거나, 계산 domain을 분할해 `REDUCE_SCATTER`하는 방법을
검토할 가치가 있다.

### 4순위: opacity fill의 memory order 개선

`kap_ion(nbin,nleaf)`은 Fortran에서 첫 index가 contiguous인데
`gas_opacity_fill()`은 `leaf`를 inner loop로 사용해 stride-`nbin` write를 한다.
cross section 배열을 사전 계산한 뒤 `leaf` outer, `bin` inner로 바꾸면 SIMD와
cache write 효율이 좋아진다.

### 5순위: build와 profile 기반 최적화

- 기본 target이 항상 `clean`을 수행하고 module dependency가 없어 매번 전체
  rebuild한다. compiler-generated dependency 또는 명시적 `.mod` dependency를
  추가하면 개발 시간이 줄어든다.
- `-ipo`는 link time과 재현성에 영향을 주므로 실제 hotspot을 `VTune`, `gprofng`,
  또는 `-qopt-report`로 확인한 뒤 유지 여부를 결정하는 편이 낫다.
- peel pass는 converged transport를 한 번 더 실행하므로 본질적으로 비싸다.
  observer 수가 많다면 동일 ray의 observer geometry 계산과 opacity walk 비용을
  별도 profile해야 한다.

## 7. 테스트 보강 제안

현재 suite는 물리 gate가 잘 분리되어 있으나 scattered peel의 자동 수치 gate가
약하다. 다음을 우선 추가하는 것이 좋다.

1. HG 산란을 `10^6`회 반복하며 매 step `abs(norm(k)-1) < 1e-12`를 검사하는 test
2. `g={-0.9999,0,0.9999}`, `cos={-1,1}`에서 phase function finite/positive 검사
3. 같은 seed의 `fast=2`/`precise` image sum 상대차 `<=1e-12` gate
4. 비정방 `(nx,ny,nz)` density cube의 세 축 centering test
5. iteration 미수렴 시 warning/metadata 존재 여부 test
6. NaN sentinel과 `is_finite`를 ifort/ifx/gfortran release flags별로 검사

`tests/peel/check_peel.py`의 현재 정량 gate는 optically-thin direct image와
zero-scattering을 검증하므로 이번에 문제가 된 scattered HG 경로를 실행하지 않는다.

## 8. 권장 적용 순서

1. 본 감사의 산란 기저/HG/read_3D 수정과 Makefile 변경 유지
2. CI에서 `peel_dbg` 또는 더 작은 scattered case를 fast/precise 양쪽으로 실행
3. convergence 종료 상태와 입력 validation 추가
4. cross section 사전 계산
5. cell solve 분할과 tally collective 개선
6. 전체 physics gate를 ifx `fast=2`, ifx `precise`, gfortran `-O3`로 비교

Intel 문서상 `precise`는 value-unsafe 최적화를 막고, `fast=2`는 더 공격적인
재결합과 contraction을 허용한다
([Intel `fp-model` reference](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2025-0/fp-model-fp.html)).
따라서 production 기본을 `fast=2`로 유지한다면 compiler version 변경 때마다
수치 regression baseline을 갱신하지 말고, 기존 tolerance를 통과하는지 먼저
확인해야 한다.

## 9. Cloudy c23.01 atomic data 비교

### 9.1 비교 범위와 해석상의 주의점

로컬 `~/CLOUDY/c23.01`의 `data/`와 실제 이를 읽는 `source/`를 함께 확인했다.
비교 대상은 photoionization, electron-impact ionization, RR/DR, H charge
transfer, recombination line, collision strength/A-value이다. rate 비교는 H II
영역의 대표값인 `T=10,000 K`에서 수행했으며, RR/DR은 `5,000`과 `20,000 K`도
추가 확인했다.

Cloudy를 단일 atomic database로 간주하면 안 된다. c23.01은 내부 H/He model,
Verner/Badnell 계수, CHIANTI, Stout, 그리고 일부 이온 전용 구현을 조합한다.
따라서 아래의 "Cloudy 값"은 `data` 파일만 비교한 것이 아니라 가능한 경우 실제
Cloudy source가 선택하는 분기까지 반영한 값이다. Cloudy 실행 결과 전체를
MoCHII와 benchmark한 것은 아니므로 thermal/ionization 구조의 최종 차이는 별도
통합 테스트가 필요하다.

### 9.2 결론과 갱신 우선순위

| 우선순위 | 항목 | 판정 |
|---|---|---|
| P0 | 고차 H charge transfer | C/N/O/Ne/S/Ar와 Si 일부에서 MOCASSIN의 한 단계 높은 경계 계수가 배정됨. 데이터 갱신 전 reaction-stage mapping 수정 필요 |
| P1 | Si II -> Si I RR+DR | Cloudy의 2023 Badnell 표가 `10^4 K`에서 현재 MoCHII의 0.476배. 원자료 확인 후 갱신 권장 |
| P1 | He I recombination lines | MoCHII 출력에 없음. H II 영역 진단용 Porter He I emissivity 추가 권장 |
| P2 | collisional ionization | Cloudy 기본은 Voronov에 Dere(2007) 비율을 곱한 hybrid. 선택 옵션으로 추가 권장 |
| P2 | Cl/Ca charge transfer | MoCHII에는 비어 있으나 Cloudy에는 일부 고차 rate가 있음. 원문 확인 후 보강 |
| P2 | Fe II/Fe III line model | MoCHII 13/14 level과 Cloudy Stout 80/70 level의 line benchmark 필요 |
| 유지 | photoionization cross section | 양쪽 모두 Verner 계열 ground-state fit. 교체 근거 없음 |
| 보완 | H I/He II recombination lines | MoCHII는 SH95 case B를 직접 사용하지만 Cloudy 주 선은 H-like model atom 결과. case 선택, 범위 검사, collisional/transfer 효과 보강 필요 |
| 유지 | CHIANTI line data | MoCHII 11.0.2가 Cloudy bundle 10.0.1보다 최신. Cloudy 파일로 downgrade 금지 |

### 9.3 Photoionization cross section

MoCHII의 `PHOTO`는 VFKY96 outer-shell fit, 자료가 없는 Cl 등은 VY95 subshell
fit을 사용한다. Cloudy의 `source/atmdat_adfa.cpp`도 같은 Verner & Yakovlev
(1995), Verner et al. (1996) `PH1/PH2` 식을 사용한다. threshold 근처 resonance나
excited-state photoionization을 상세하게 다루려면 Opacity Project/level-resolved
자료가 필요하지만, 현재 MoCHII의 ground-state continuum 모델을 Cloudy c23.01
계수로 바꿔 얻는 이점은 없다.

판정: **현행 유지**. 향후 He I excited-state 또는 metastable population을 명시적으로
풀 때만 level-resolved cross section을 별도 도입한다.

### 9.4 Radiative/dielectronic recombination

Cloudy의 `badnell_rr.dat`와 `badnell_dr.dat`는 각각 `20230511`, `20230512`
magic date를 갖는다. MoCHII는 CHIANTI 11.0.2의 `.rrparams/.drparams`를 변환한다.
공통으로 Badnell 형식을 사용하는 이온의 총 `RR+DR`을 직접 평가한 결과, C, N,
O, Ne, Mg, 대부분의 S/Ar/Si/Cl/Ca 이온은 파일 유효숫자 범위에서 동일했다.

의미 있는 예외는 `Si+ + e -> Si0`이다.

| T [K] | MoCHII [cm3 s-1] | Cloudy Badnell [cm3 s-1] | Cloudy/MoCHII |
|---:|---:|---:|---:|
| 5,000 | 4.038e-12 | 2.678e-12 | 0.663 |
| 10,000 | 3.416e-12 | 1.626e-12 | 0.476 |
| 20,000 | 5.558e-12 | 4.708e-12 | 0.847 |

CHIANTI의 Si II 파일은 2012 Abdel-Naby et al. 계수이고, Cloudy c23.01의 2023
Badnell 집계표는 일부 fit coefficient가 바뀌었다. Si I/Si II ion fraction과
`[Si II]` line은 이 차이에 민감할 수 있으므로 P1 갱신 대상으로 둔다. Cl V ->
Cl IV는 `10^4 K`에서 2.5% 차이, S II -> S I은 0.3% 이하로 우선순위가 낮다.

Cloudy Badnell 표에 없는 저이온 Ar, Ca, Fe에 대해 Cloudy는 RR에 Verner fallback,
DR에 isoelectronic mean 또는 전용 fit을 사용한다. 비교 계산에서 Badnell 항이 0인
것을 Cloudy의 총 recombination이 0이라는 뜻으로 해석해서는 안 된다. 이 이온들은
Cloudy의 `save recombination coefficients` 출력과 MoCHII `alpha_rec()`를 직접
비교하는 후속 benchmark가 필요하다.

### 9.5 Collisional ionization

MoCHII는 Voronov (1997)를 그대로 사용한다. Cloudy c23.01의 실제 default는
`source/init_defaults_preparse.cpp`에서 `HYBRID`이며, Voronov rate에
`source/atmdat_adfa.cpp`의 Dere (2007)/Voronov 비율을 곱한다. MoCHII가 추적하는
stage 중 큰 보정은 다음과 같다.

| ionization | Cloudy hybrid/Voronov |
|---|---:|
| Mg I -> Mg II | 0.3793 |
| Si III -> Si IV | 0.4492 |
| Cl I -> Cl II | 0.5412 |
| Si I -> Si II | 0.7328 |
| S I -> S II | 1.3572 |
| Fe III -> Fe IV | 1.8240 |

보정 폭은 크지만 `10^4 K`의 정상 H II 영역에서는 금속 ionization threshold의
Boltzmann suppression 때문에 photoionization과 charge transfer가 대개 우세하다.
반면 저이온화 전선, PDR, shock/hot cell에서는 차이가 커질 수 있다.

판정: 기존 Voronov를 삭제하지 말고 `ci_model = voronov|dere_hybrid`와 같은 선택
옵션으로 구현한다. H II benchmark가 안정되면 hybrid를 기본값으로 바꿀 수 있다.

### 9.6 Charge transfer: EXHALE 원본과 MoCHII 이식 오류

#### 9.6.1 Huang 원문과 EXHALE

로컬 원문
`~/RT_Codes/Exoplanetary_Atmospheres/references/Huang_2023_ApJ_951_123.pdf`의
Table 4(pp. 25-26)를 렌더링해 계수, `T4`의 곱/나눗셈 위치와 반응물을 EXHALE
`src/modules/radiation/charge_exchange.f90`과 직접 대조했다. 정적 대조에서 다음은
올바르다.

- EXHALE의 23개 metal-H/H+ 기본 반응과 40개 He-H/metal-He/metal-metal 선택 반응은
  Table 4에 실제로 인쇄된 63개 행을 모두 포함한다.
- `kf`, `gp`, `p4`가 각각 대괄호 안의 `exp(-a*T4)`, activation factor
  `exp(-E/T4)`, `ln(T)` polynomial을 올바르게 구분한다. 63개 식에서 계수나 지수
  부호의 불일치를 찾지 못했다.
- donor는 `stage s -> s+1`, acceptor는 `stage s -> s-1`로 명시돼 있다. 예를 들어
  A3 `Mg+ + H+`와 A4 `Mg2+ + H`는 모두 Mg+/Mg2+ 경계에 반대 방향으로 들어가며,
  residual/Jacobian의 부호도 이 정의와 일치한다.
- 이 구현과 반응 descriptor는 2026-06-11의 EXHALE commit `10c248a`부터 있었고,
  이후 변경은 주로 OpenMP thread-local 상태 처리였다. MoCHII 최초 import
  `96a4066`보다 먼저 존재했다.

논문 본문 Section 2.4.2는 "65 reactions"라고 쓰지만, 실제 Table 4의 reactant 행은
63개다. EXHALE의 `23 + 2 + 6 + 32 = 63`은 **인쇄된 표에는 완전**하다. 빠진 것으로
보이는 두 반응을 원문만으로 특정할 수 없으므로, EXHALE 문서가 그 차이를
"Groups C-D의 두 행"이라고 단정한 부분은 근거가 부족하다. 이는 발견된 rate-code
오류라기보다 Huang 논문의 본문/표 불일치이며, 저자 또는 supplement 확인 전에는
임의의 두 반응을 추가하면 안 된다.

EXHALE에 남은 위험은 반응별 golden/unit test가 없다는 점이다. 광범위한 모델
regression은 있으나 `cx_rate(1:63,T)`와 stage/source-term 부호를 독립적으로 고정한
테스트는 찾지 못했다. 또한 각 원자료의 유효 온도 범위를 반응별로 보존하지 않고
최종 rate만 `[0,1e-6]`으로 자르므로, 모델 온도 범위가 넓어질 때 fit 외삽을 경고하지
못한다. 이 두 항목은 보강 대상이지만, 이번 대조에서 EXHALE의 Table-4 전사 또는
stage 배정 오류는 확인되지 않았다.

#### 9.6.2 MoCHII가 정확히 가져온 부분

MoCHII `tools/fitting/make_element_data.py:64-89`의 `CX`는 첫 ionization 경계에
들어가는 H/H+ 반응이다. C(A15/A16), N(A17/A18), O(A13/A14), S(A19/A20),
Mg(A1/A2), Fe(A5/A6)는 EXHALE/Huang 식과 동일하며, `species_mod.f90`도
`CXI*n(H+)`를 ionization, `CXR*n(H0)`를 recombination에 넣으므로 방향이 맞다.
따라서 이 부분은 EXHALE에서 가져오면서 망가진 것이 아니다.

다만 평가 범위는 같지 않다. MoCHII의 KF96 form 6/7은 `T4`를 `0.5-5`로 조용히
clamp하고, EXHALE은 원식을 현재 cell 온도에서 평가한 뒤 rate만 `[0,1e-6]`으로
제한한다. 따라서 공통 계수도 `5000-50000 K` 밖에서는 서로 다른 값이 된다. 이는
계수 전사 오류는 아니지만, 두 코드 비교 시 반드시 기록해야 할 extrapolation
정책 차이다.

Mg와 Fe의 두 번째 경계에서 recombination 방향 A4 `Mg2+ + H`와 A8
`Fe2+ + H`도 올바른 transition 2에 있다. 그러나 반대 방향 A3 `Mg+ + H+`와 A7
`Fe+ + H+`는 MoCHII에 빠졌다. `10^4 K`에서 각각 `1.393e-14`, `1.145e-10
cm3 s-1`이므로, 특히 Fe+/Fe2+ 경계에는 무시 가능하다고 미리 가정하면 안 된다.

#### 9.6.3 MOCASSIN 계수를 섞으며 생긴 stage 오류

MoCHII의 `CX_HIGH`는 EXHALE이 아니라
`~/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/source/update_mod.f90:1395-1445`의 오래된
Kingdon & Ferland 계수에서 왔다. MOCASSIN은 `ionRatio(element,ion) =
X(ion+1)/X(ion)`을 만들 때 같은 `chex(element,ion)`을 recombination denominator에
넣는다. 따라서 `chex(C,2)`는 C+/C2+ 경계, `chex(C,3)`은 C2+/C3+ 경계다. 주석의
`C+2`는 생성되는 lower ion을 표시하며, 그 계수를 C2+ 반응물의 C2+ -> C+ 반응으로
읽으면 한 단계 낮은 경계에 배정된다.

현재 MoCHII는 C/N/O/Ne/S/Ar에서 transition 2에 MOCASSIN `chex(*,3)`을,
C와 S의 transition 3에는 `chex(*,4)`를 넣었다. Si도 `chex(14,2:4)`를 각각
transition 1:3에 넣어 모두 한 경계 아래로 밀렸다. 반면 Mg/Fe는 Huang label을
별도로 적용해 transition 번호가 맞다. 즉 `CX_HIGH` 전체를 일괄 `+1` 또는 `-1`
shift해서는 안 되고, 반응별 source/product charge를 다시 써야 한다.

대표적인 `10^4 K` rate는 다음과 같다. 화살표는 H0와의 recombination 방향이다.

| reaction | 현재 MoCHII | 올바른 비교값 | 판정 |
|---|---:|---:|---|
| C2+ -> C+ | 3.273e-9 | Cloudy 1.036e-12 | MOCASSIN의 C3+ -> C2+ 계수를 한 단계 낮게 사용 |
| N2+ -> N+ | 3.328e-9 | Cloudy 9.703e-10 | 한 단계 높은 MOCASSIN row 사용 |
| O2+ -> O+ | 4.142e-9 | Cloudy 8.598e-10 | 한 단계 높은 MOCASSIN row 사용 |
| Ne2+ -> Ne+ | 5.655e-9 | Cloudy 1.0e-14 placeholder | 한 단계 높은 MOCASSIN row 사용 |
| S2+ -> S+ | 2.298e-9 | Cloudy 1.0e-14 placeholder | 한 단계 높은 MOCASSIN row 사용 |
| Ar2+ -> Ar+ | 4.399e-9 | Cloudy 1.0e-14 placeholder | 한 단계 높은 MOCASSIN row 사용 |
| Si+ -> Si0 | 5.145e-9 | Huang A10 2.371e-12 | 실제로는 A12와 유사한 Si2+ -> Si+ 식을 사용 |
| Si2+ -> Si+ | 4.099e-10 | Huang A12 5.238e-9 | MOCASSIN의 Si3+ -> Si2+ 식을 사용 |
| Si3+ -> Si2+ | 7.714e-9 | Cloudy 4.099e-10 | MOCASSIN의 Si4+ -> Si3+ 식을 사용 |

Cloudy의 `1e-14` 값은 정밀 측정값이 아니라 "매우 작음" placeholder에 가깝다.
Cloudy와 독립적인 확인으로 MOCASSIN의 실제 ratio loop와 CMacIonize의
`ION_C_p2` mapping도 같은 경계 해석을 사용한다. 따라서 위 판정의 핵심은 Cloudy
rate를 정답으로 가정한 것이 아니라, 원본 코드에서 계수가 사용되는 charge boundary를
추적한 결과다.

Si는 MoCHII 추가 commit `e106203`에서 별도로 들어왔으며 EXHALE A9-A12를 제대로
가져오지 않았다. transition 1에는 A9 `Si + H+`/A10 `Si+ + H`, transition 2에는
A11 `Si+ + H+`/A12 `Si2+ + H`가 들어가야 한다. 현재는 forward가 모두 0이고
MOCASSIN recombination 식만 한 단계씩 낮게 배정돼 있다. 이는 EXHALE 원본 문제가
아니라 MoCHII 이식 과정의 오류다.

#### 9.6.4 수정 순서

1. 각 CX row를 `(donor element, donor charge, acceptor element, acceptor charge,
   source)` 구조로 생성하고 `TRANSITION`은 여기서 계산한다. 배열 위치나 ion 주석을
   반응물로 해석하지 않는다.
2. MoCHII가 추적하는 원소에 대해 Huang A1-A20을 기준으로 다시 생성한다. 기존
   C/N/O/S/Mg/Fe 첫 경계는 유지하고, Mg A3, Fe A7과 Si A9-A12를 보강한다.
3. Huang 표에 없는 고차 반응은 KF96 원표 또는 Cloudy c23.01의 실제 선택값에서
   source/product charge를 명시해 다시 생성한다. 특히 Cloudy는 N2+/O2+에 2006년
   이후 전용 fit을 덮어쓴다.
4. `T={5e3,1e4,5e4 K}`에서 reaction별 golden table을 만들고 generator 출력과
   Fortran `cx_rate()`를 함께 검사한다. EXHALE에도 같은 63-row test를 추가한다.
5. 수정 전후 C/N/O/Ne/S/Ar/Si/Mg/Fe ion fraction과 주요 line luminosity를 같은
   Cloudy H II benchmark에서 비교한다.

이 항목은 ion fraction을 수 배에서 수천 배까지 바꿀 수 있으므로 **P0**이다.
본 감사에서는 오류를 문서화했지만 atomic 파일은 아직 변경하지 않았다. 올바른
원자료 선택과 Cloudy 통합 benchmark를 먼저 고정해야 physics regression을 해석할
수 있다.

### 9.7 Collision strength, A-value, line model

MoCHII 생성원은 CHIANTI `11.0.2`, Cloudy c23.01 bundle은 `10.0.1`이다. 따라서
MoCHII의 CHIANTI 기반 `nlevel_*`와 `cooling_tier1_*`를 Cloudy bundle로 교체하면
downgrade가 된다. 현재 1e3-1e5 K Chebyshev fit 오차도 대부분 2% 이내이므로 fit
정밀도 때문에 갱신할 항목은 보이지 않았다.

다만 Cloudy는 저이온 H II 진단종에서 CHIANTI만 사용하지 않는다. 기본 masterlist
기준으로 C I-III, O I-III, Ne I-III, S I-III, Ar I-IV, Si I-IV, Fe I-III 등은
주로 Stout를 사용하고, N II/III, Mg II, Cl II-IV, Ca II/V 등은 CHIANTI를 선택한다.
따라서 "CHIANTI 11이 더 최신"이라는 사실만으로 line emissivity가 Cloudy보다 더
정확하다고 결론내릴 수 없다.

특히 MoCHII Tier-2는 Fe II 13 level, Fe III 14 level인 반면 Cloudy Stout default는
최소 80/70 level이며 Fe II에는 Bautista (2018), Tayal (2018), Smyth (2019) 계열
자료가 포함된다. MoCHII Tier-1 총 cooling fit은 더 많은 CHIANTI excitation을
압축하지만, 개별 Fe line spectrum과 density coupling은 이 차이를 흡수하지 못한다.

판정: 주요 optical/IR 진단선 `[O III] 4363/5007`, `[N II] 5755/6584`,
`[S II] 6716/6731`, `[Ne III] 15.5 um`, `[Ar III/IV]`, `[Fe II/III]`를 동일한
`T, ne, ion density`에서 Cloudy와 비교한다. 차이가 큰 이온만 level 수 또는 원자료를
선택적으로 보강하고, 전체 CHIANTI 자료를 교체하지 않는다.

### 9.8 H/He recombination line data

MoCHII의 H I와 He II 출력은 Storey & Hummer (1995; SH95) case B 표를 직접
보간한다. Cloudy도 같은 표를 `Ca A`/`Ca B` reference line으로 제공하지만, 일반
출력의 `H  1`/`He 2` 선은 level-resolved H-like model atom으로 다시 계산한다.
따라서 원자료가 같다는 이유만으로 실제 예측선이 같다고 볼 수 없다. H I에 대한
상세 호출 경로와 수치 비교는 10절에 정리했다.

Cloudy는 이와 별도로 `he1_case_b.dat`와 He-like model을 통해 Porter et al.
(2012, 2013)의 He I emissivity를 제공한다. MoCHII는 He ionization balance,
free-bound/free-free continuum, He I two-photon은 계산하지만 He I recombination
line map은 출력하지 않는다.

H II 영역 연구 코드라면 He I 4471, 5876, 6678, 7065, 10830 A는 helium abundance,
temperature, density, optical-depth 진단에 직접 필요하다. 따라서 Porter grid를
읽는 보간 모듈과 line output을 P1 기능으로 추가하되, 단순 case-B table만 넣을지
collisional/metastable 및 optical-depth correction까지 풀지는 별도 옵션으로
구분해야 한다.

### 9.9 권장 실행 순서

1. CX reaction-stage golden table 작성 및 P0 stage 오류 수정
2. 표준 H II sphere를 Cloudy/MoCHII에서 동일 abundance/SED/density로 실행해 ion
   fraction과 주요 line ratio 비교
3. Cloudy 2023 Badnell Si II RR/DR을 별도 branch에서 적용하고 Si line 영향 측정
4. Porter He I case-B emissivity 출력 추가
5. Dere hybrid collisional-ionization 선택 옵션 추가
6. Fe II/III level 확장 여부를 line benchmark 결과로 결정

atomic data 파일에는 source version, 원 논문, 생성 날짜뿐 아니라 source checksum과
각 reaction의 명시적 charge도 기록해야 한다. 현재처럼 배열 위치와 주석에 stage
의미를 맡기면 같은 종류의 전사 오류를 자동으로 잡기 어렵다.

## 10. Cloudy와 MoCHII의 수소 recombination line 계산

### 10.1 결론

Cloudy c23.01의 일반 Balmer line은 SH95 case-B emissivity를 그대로 보간한 값이
아니다. SH95와 독립적인 level-resolved recombination coefficient를 각 H 원자
준위의 source term으로 넣고, radiative cascade, electron/heavy-particle collision,
continuum pumping, excited-state photoionization, induced/three-body recombination,
line escape와 destruction을 포함한 population matrix를 푼 뒤 선을 계산한다.

반면 MoCHII는 수렴한 각 cell에서

```text
epsilon_ul = n_e n(H+) E_SH95_caseB(T, n_e)  [erg s^-1 cm^-3]
```

를 적용하는 post-processing이다. SH95 표 자체의 density-dependent cascade와
collisional l-mixing은 포함하지만, 그 cell의 H I population, radiation field,
Balmer optical depth와 line destruction은 emissivity에 되먹임되지 않는다.

따라서 전형적인 저밀도, Lyman-thick H II 영역에서는 두 값이 수 % 이내로 가까울
수 있으나, partially ionized front, 높은 온도/밀도, 강한 continuum pumping 또는
Balmer line trapping이 있는 경우 MoCHII의 H line은 Cloudy의 `H  1`보다 단순한
case-B reference이다. 동일한 물리량을 비교하려면 먼저 MoCHII를 Cloudy의 `Ca B`와
비교하고, 실제 model 차이는 `H  1`과 별도로 비교해야 한다.

### 10.2 Cloudy에서 실제로 계산되는 경로

Cloudy source에서 확인한 경로는 다음과 같다.

1. `iso_create.cpp`가 `iso_recomb_setup(ipH_LIKE)`를 호출해
   `data/h_iso_recomb.dat`의 `(n,l)`-resolved radiative recombination coefficient를
   읽는다. 이것은 line emissivity 표인 `HS_e1b.dat`와 다른 입력이다.
2. `iso_radiative_recomb.cpp`는 resolved level에는 이 coefficient를 온도 보간해
   넣고, collapsed level에는 hydrogenic coefficient를 사용한다. 유한한 model
   atom에서 빠진 높은 n의 recombination은 최상위 collapsed level에 `topoff`한다.
   기본 H atom은 n=10까지 l-resolved이고, 그 위 15개 n-shell은 collapsed level이다
   (`iso.cpp:190-197`).
3. 같은 routine은 ground-level을 제외한 합을 `RadRec_caseB`로 따로 저장하지만,
   일반 level source에는 각 continuum의 net escape probability를 곱한다.
   `case B` command를 지정한 경우에만 ground recombination의 net escape를 사실상
   0으로 둔다.
4. `iso_ionize_recombine.cpp`의 continuum-to-level rate에는 radiative/dielectronic/
   induced/three-body recombination이 들어간다. 반대 방향에는 photoionization,
   collisional ionization과 secondary ionization이 들어간다.
5. `iso_level.cpp`가 `RateCont2Level*n(H+)`를 source vector로 사용하고, 모든 level
   사이의 collisional excitation/de-excitation, `A_ul*P_loss`, continuum pumping과
   관련 sink/source를 포함한 행렬을 푼다. 따라서 neutral H의 ground-state
   collisional excitation도 Balmer population에 기여할 수 있다.
6. `prt_lines_hydro.cpp`가 l-resolved component를 같은 n-to-n' transition으로
   합친다. `lines_service.cpp:set_xIntensity()`의 gross local line emissivity는
   `n_u A_ul P_esc h nu`이며, continuum subtraction을 요청하면 stimulated/background
   항까지 반영한 observed emissivity를 별도로 만든다.

Cloudy의 `case B` command도 곧바로 SH95 표 lookup으로 계산 방식을 바꾸는 명령은
아니다. `parse_caseb.cpp`에 따르면 ground recombination과 Lyman optical depth를
case-B 쪽으로 강제할 뿐 model atom solve는 계속한다. 엄격한 SH95 재현 시험은
`case B hummer no photoionization no pdest`, `no induced processes`, 많은 resolved/
collapsed level 등의 조건을 함께 사용한다. `limit_caseb_h_hs87.in`은 이 조건에서
주 `H  1`과 `Ca B` H-beta를 2% 이내로, `limit_caseb_h_den4_temp4.in`은 여러 H line을
SH95 대비 5% 이내로 검사한다.

### 10.3 Cloudy line label의 의미

Cloudy 출력에서 같은 H-beta 파장 부근에 여러 line이 존재한다.

| label | 의미 | MoCHII와의 대응 |
|---|---|---|
| `H  1` 4861.32 A | full H-like model atom의 주 예측값 | 직접 대응 없음 |
| `Ca B` 4861.32 A | `HS_e1b.dat`의 SH95 case-B 직접 보간 | 현재 MoCHII와 대응 |
| `Ca A` 4861.32 A | SH95 case-A 직접 보간 | MoCHII에 없음 |
| `CaBo` 4861.32 A | Ferland (1980)의 오래된 analytic case-B 근사 | 사용하지 말 것 |
| `Q(H)` 4861.32 A | ionizing photon rate와 covering factor로 추정한 H-beta | global sanity check |

Cloudy 기본 line list는 2000 A보다 긴 파장에 air wavelength를 쓰므로 H-alpha/H-beta가
6562.80/4861.32 A로 표시된다. MoCHII의 6564.6/4862.7 A는 vacuum wavelength이다.
이는 atomic energy 차이가 아니므로 transition `(n_u,n_l)` 또는 같은 wavelength
convention으로 비교해야 한다.

### 10.4 MoCHII 구현과 내부 불일치

`lines_mod.f90:63-112`는 각 leaf의 `n_e`, `n_H(1-x_HI)`와
`sh95_emis()`를 곱해 H-alpha, H-beta, H-gamma, H-delta, Pa-alpha, Pa-beta,
Br-gamma luminosity와 emissivity cube를 만든다. `sh95_mod.f90`은 log T-log n_e에서
log emissivity를 bilinear interpolation한다.

중요한 문제는 H ionization balance의 case 선택과 line case가 독립이라는 점이다.
`ion_balance_mod.f90`은 `par%case_ab`에 따라 Hui & Gnedin의 alpha_A 또는 alpha_B를
사용하고, explicit diffuse field는 setup에서 case A를 강제한다. 그러나
`lines_mod.f90`은 어떤 설정에서도 SH95 **case B**만 사용한다. 따라서
`diffuse_field=.true.`/`case_ab='A'` 모델도 line output은 case B이며, metadata에는
이 물리적 혼합이 충분히 드러나지 않는다.

또한 SH95 H I grid는 `T=300-30000 K`, `n_e=10^2-10^14 cm^-3`인데,
`sh95_emis()`는 양쪽 범위를 벗어난 값을 모두 경계로 조용히 clamp한다. MoCHII의
기본 thermal upper bound는 50000 K이므로 30000-50000 K cell에서 H emissivity가
30000 K 값으로 고정될 수 있다. Cloudy의 SH95 helper는 낮은 density만 최저점으로
clamp하고, 온도 범위 밖 또는 최고 density 초과 시 `-1`을 반환해 해당 reference가
유효하지 않음을 표시한다.

### 10.5 동일 SH95 표의 수치 대조

MoCHII의 `sh95_hi_caseB.txt`와 Cloudy의 `HS_e1b.dat`를 공통 grid 10x13의 모든
점에서 직접 비교했다. 두 파일은 같은 SH95 자료의 서로 다른 repackaging이다.
7개 line에서 최대 상대차는 0.11-0.17%였고, 이는 주로 저장 자릿수 차이다.

| T [K] | n_e [cm^-3] | line | MoCHII | Cloudy `Ca B` | 상대차 |
|---:|---:|---|---:|---:|---:|
| 500 | 1e2 | H-alpha | 4.545e-24 | 4.542e-24 | +0.066% |
| 5000 | 1e4 | H-beta | 2.224e-25 | 2.224e-25 | 0.000% |
| 10000 | 1e4 | H-alpha | 3.531e-25 | 3.530e-25 | +0.028% |
| 10000 | 1e4 | H-beta | 1.240e-25 | 1.240e-25 | 0.000% |
| 20000 | 1e4 | H-beta | 6.589e-26 | 6.589e-26 | 0.000% |
| 30000 | 1e14 | Br-gamma | 1.212e-27 | 1.211e-27 | +0.083% |

계수 단위는 모두 `4 pi j/(n_e n_p) [erg cm^3 s^-1]`이다. 예를 들어
`T=10000 K`, `n_e=1e4 cm^-3`에서 MoCHII의 case-B H-alpha/H-beta는 2.8476이고
Cloudy `Ca B`는 2.8468이다. 즉 현재 표를 Cloudy 표로 교체할 필요는 없다.

그러나 Cloudy regression suite 자체에서도 full model과 `Ca B`는 항상 같지 않다.
순수 case-B 시험에서는 H-beta가 약 0.6% 수준에서 일치하지만, `nlr_paris.in`과
`nlr_lex00.in`에서는 `Ca B/H 1`이 각각 0.921, 0.929로 full H-beta가 case-B
reference보다 약 8.6%, 7.6% 크다. 반대로 optically thick continuum-emission
시험 `h_t4_conemis_thick.in`에서는 `Ca B/H 1=1.075`로 full H-beta가 약 7.0%
작다. 이 수치는 보편적 correction factor가 아니라, 빠진 물리의 효과가 모델에
따라 부호도 바뀐다는 예시다.

### 10.6 권장 수정과 검증

1. **P0: 범위 밖 clamp를 중단한다.** 최소한 T/ne가 표 밖이면 warning과 output
   metadata를 남기고, `T>30000 K`는 H line 계산에서 제외하거나 명시적 extrapolation
   정책을 선택한다. H II 연구의 일반 범위에서는 무음 clamp보다 fail-fast가 낫다.
2. **P1: line case를 명시한다.** SH95 case-A 표도 읽어 `par%case_ab`와 맞추거나,
   `h_line_case='A'/'B'`를 독립 옵션으로 두되 ionization case와 다를 때 경고한다.
   output header/HDF5 metadata에도 `SH95 reference, case B, intrinsic`을 기록한다.
3. **P1: collisional H I component를 추가한다.** 우선 H-alpha/H-beta에 대해
   ground 및 n=2에서의 excitation을 별도 component로 출력하면 ionization front와
   고온 cell에서의 bias를 정량화할 수 있다. recombination과 합치기 전에 Cloudy의
   `H  1 - Ca B` 및 local H0-weighted emissivity와 benchmark한다.
4. **P2: H model atom/escape 옵션을 검토한다.** Balmer optical depth, pumping 또는
   고밀도 영역이 연구 대상이면 fixed correction보다 작은 H n-level solver와 escape
   probability를 도입해야 한다. 전형적 저밀도 H II 영역만 목표라면 이 복잡도는
   기본값이 아니라 선택 기능으로 두는 것이 적절하다.
5. **회귀 테스트를 둘로 나눈다.** `(T,n_e)` 단위 테스트는 MoCHII 대 Cloudy
   `Ca A/Ca B`를 0.2% tolerance로 검사하고, 통합 H II sphere는 Cloudy `H  1`과
   H-alpha/H-beta, H-beta luminosity 및 front별 emissivity를 비교한다. 두 비교를
   섞으면 table transcription 오류와 model-physics 차이를 구분할 수 없다.

이번 수치 비교는 Cloudy executable을 새로 빌드해 실행한 것이 아니라, c23.01의
실제 source/data 호출 경로와 배포 regression input의 monitor 기준값을 조사하고
두 SH95 data file을 직접 파싱한 결과이다. full spatial benchmark는 위 5번의 별도
검증 항목으로 남는다.
