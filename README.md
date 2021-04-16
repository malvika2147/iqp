# iqp

Simulates the IQP protocol for g=0
- The verifier generates a P_s such that P_s^TP_s = 0 and then adds many random rows orthogonal to s, to get P.
- The adversarial prover sends samples of the form y(u,v) = \sum_{p \in P} p <u,p><v,p> and these samples are always orthogonal to s.
- This code checks the above on random P (generated such that P_s^TP_s = 0) and checks how many prover samples are unique and non-zero.
