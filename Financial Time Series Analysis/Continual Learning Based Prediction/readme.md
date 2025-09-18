# Continual Learning for Financial Time Series Forecasting

## Introduction

Continual learning addresses the challenge of adapting deep learning models to evolving financial market data without catastrophic forgetting of prior knowledge. This is critical in financial time series forecasting where data distribution and market regimes shift continuously.

## Methods Implemented

- **Elastic Weight Consolidation (EWC)**  
  Uses a regularization term weighted by the Fisher Information Matrix to constrain important parameters from changing drastically, enabling retention of past knowledge while learning new tasks.

- **Learn Without Forgetting (LWF)**  
  Uses knowledge distillation loss by encouraging the model to mimic outputs from a previously trained version (teacher), combining this with ground truth loss for balanced adaptation.

- **GDumb**  
  A replay-based approach that maintains a balanced memory buffer of past samples and retrains from scratch on stored data, without assumptions about task boundaries or label distribution shifts.

## Baseline Approaches

- **Naive Training**: Sequential training without any forgetting mitigation.
- **Joint Training**: Retraining on the combined old and new data.
- **Re-incremental Training**: Training only on newly added data without past knowledge.

## Evaluation

The methods are evaluated on Bitcoin price datasets with incremental training windows. 

## References

- Kirkpatrick et al., 2017 (EWC)  
- Li and Hoiem, 2017 (LWF)  
- Prabhu et al., 2020 (GDumb)  

---

