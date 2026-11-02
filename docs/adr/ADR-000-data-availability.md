# ADR-000 데이터 가용성 확인

- 날짜: 2026-11-02
- 상태: 확정
- 관련 이슈: #1

## 문제

KRAS G12D 결합 예측 모델(GCN)을 학습시킬 데이터가 충분한지,
scaffold split 평가가 유효한지, assay 조건을 구분할 수 있는지
확인이 필요했다. 이 결과에 따라 전이학습 도입 여부와
룩업 백업 사용 여부가 갈린다.

## 데이터 출처

BindingDB 202611 (2026-10-31 release)
- `BindingDB_All.tsv` (8.4 GB, 3,243,660행), md5 검증 완료
- `BindingDB_Assays.tsv` (assay 설명)
- `BindingDB_rsid_eaids.tsv` (조인 키)

## 필터링

| 단계 | 건수 |
|---|---|
| 전체 | 3,243,660 |
| Target Name에 KRAS 계열 포함 | 30,278 |
| 변이 G12D | 8,284 |
| 복합체 제외 (단독 KRAS만) | 4,956 |
| IC50 측정, 관계기호 `=` | 3,934 |
| 이상치 제외 (pIC50 > 10) | **3,881** |

복합체 제외 기준: Target Name에 SOS1, RAF, PI3K 등 파트너
단백질이 함께 표기된 것. 단백질-단백질 상호작용 억제를 측정한
것으로 기전이 다름.

관계기호: `=` 24,430건만 회귀에 사용. `>` 2,620건, `<` 1,529건은
별도 보관(분류용으로 활용 가능).

## 데이터 특성

### 구조 다양성
- 고유 분자: 2,271
- 고유 Murcko scaffold: 938
- 단일 분자 scaffold: 651 (69%)
- 상위 10개 scaffold가 분자의 24.3% 커버
- 2개 이상 특허에 등장하는 scaffold: 56 / 938 (6%)

### 활성 분포
- pIC50 중앙값 7.20, 표준편차 1.45, 범위 3.72~10.0
- 이봉 분포 (5.1~6.5 구간과 7.9~9.2 구간에 집중)
  → 특허별 화합물 시리즈 차이로 추정

### 출처 구성
- G12D 단독 데이터는 100% US Patent (ChEMBL/논문 유래 0건)
- 고유 특허: 24개
- **상위 10개 특허가 분자의 90.7%를 차지**
- 특허당 scaffold 수는 다양 (최대 225개, scaffold 1개뿐인
  특허는 24개 중 1개)

### assay 기전 분류
3,881건 중 3,359건(87%)에 assay 설명이 연결됨.

| 기전 | 대표 assay 명 | 비고 |
|---|---|---|
| 뉴클레오타이드 교환 억제 | Coupled Nucleotide Exchange, TR-FRET | SOS1 매개 |
| 단순 결합 | KRas G12D Binding Assay | |
| 세포 내 ERK 인산화 억제 | Inhibition of KRas G12D-mediated Phosphorylation of ERK | 444건, **세포 assay** |
| 분류 모호 | Biochemical Assay | 설명 개별 확인 필요 |

→ "생화학적 IC50"이 단일 기전이 아님이 확인됨. 세부 분류는
   후속 작업으로 진행.

## 이상치 처리

pIC50 > 10 (IC50 < 0.1 nM) 6건 제외.

근거:
- 생화학 assay의 단백질 농도(확인된 사례에서 5 nM)보다 낮은
  IC50은 정확히 측정되지 않음 (tight-binding limit)
- 2개 특허(US20230279025, US20230339981)에만 집중
- US20230279025 Example 463~466은 연속 화합물인데
  0.001~0.012 nM 범위에서 값이 불규칙 → 실측보다 상한 표기로 추정

하한(pIC50 3.72, 4.30 각 1건)은 유지. 50~190 µM은 측정 가능
범위이며, 약한 화합물도 학습에 필요함.

## 결정

**원안 유지.**

1. 전이학습 도입하지 않음 — 데이터 규모 충분
2. 룩업 백업 사용하지 않음 — 실제 GCN 학습 가능
   (인터페이스는 교체 가능하게 유지)
3. **주 평가: scaffold split** — 938개 scaffold, 특허 간
   겹침 6%로 누수 위험 낮음
4. **보조 평가: 특허 split 추가** — 24개 특허를 학습/평가로
   분할. "본 적 없는 출원인의 화합물을 맞히는가"를 측정
5. 베이스라인: fingerprint + Random Forest
6. assay 기전별 분리를 후속 작업으로 진행. 세포 assay
   (ERK 인산화 444건)는 회귀 학습셋에서 제외 예정

## 한계

1. **단일 출처(특허), 동료 심사 없음.** 출원인이 자사 화합물을
   유리하게 보고할 동기가 있음
2. **출원인 집중.** 상위 10개 특허가 90.7%를 차지 →
   모델이 특정 화학 계열에 과적합될 위험. 이는 하네스의
   적용 범위(applicability domain) 검증이 필요한 직접적 근거가 됨
3. **assay 조건 미통일.** 기전이 최소 3종으로 갈림. 13%는
   assay 설명 자체가 없음
4. **논문 기반 sanity check 앵커 부재.** MRTX1133 등은 ChEMBL
   유래이며 본 데이터셋에 없음 → ChEMBL API로 별도 확보 필요
5. 특허 split의 평가셋 크기가 작을 수 있음 (특허 24개)

## 재현
data/raw/BindingDB_All.tsv (md5: 932346d34d4a3ab823b58ae6ff63b7e8)
notebooks/02_extract_kras.py → data/processed/kras_raw.csv
notebooks/03_classify.py → data/processed/kras_classified.csv
notebooks/05_patent_check.py → data/processed/g12d_unique.csv
notebooks/06_assay_join.py → data/processed/g12d_with_assay.csv



## 인용

Gilson MK, Liu T, Baitaluk M, Nicola G, Hwang L, Chong J.
BindingDB in 2015: A public database for medicinal chemistry,
computational chemistry and systems pharmacology.
Nucleic Acids Research 44:D1045-D1053. doi:10.1093/nar/gkv1072
