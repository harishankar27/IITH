# Deep Learning and Continual Learning Approaches for Financial Time Series Analysis

This repository contains implementations and experiments on forecasting financial time series data using various deep learning methods including LSTM, GRU, Transformer, and Hybrid models. The goal is accurate prediction of stock prices and financial trends by leveraging advanced neural network architectures.

Financial markets are inherently nonlinear and dynamic. Traditional linear and statistical models often fall short in capturing market complexities and adapting to shifting regimes. Deep learning methods provide capabilities for learning complex temporal dependencies and multiple features beyond historical prices.

## Models Implemented

- **LSTM**: Two-layer LSTM network capturing long-term sequential dependencies.
- **GRU**: Three-layer Gated Recurrent Unit network with dropout regularization.
- **Transformer**: Single encoder layer with multi-head attention for efficient long-term modeling.
- **Hybrid**: Combination of LSTM and GRU layers for leveraging complementary strengths.

## Continual Learning Framework

The repository also explores continual learning strategies to mitigate catastrophic forgetting when models learn sequentially from evolving market data streams:

- Elastic Weight Consolidation (EWC)
- Learn Without Forgetting (LWF)
- GDumb (Replay-based method)
- Baseline approaches: Naive training, Joint training, Re-incremental training

These methods enable incremental learning while retaining past knowledge crucial for financial predictions.

## Datasets

- Tata Motors daily stock closing prices (2010-2023) for comparative model analysis.
- Bitcoin price data sampled every 5 minutes for 60 days for continual learning evaluations.

---

