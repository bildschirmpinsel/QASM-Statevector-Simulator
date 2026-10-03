## Aim of this Project
This repository is a C++ rewrite of a quantum statevector simulator written in Python for the 2024/25 practical course `Experimental Evaluation and Characterization of Quantum Computing Systems` held at Technical University Munich.

The simulator currenly supports OpenQASM circuits with gates X, Y, Z, H, S, T, AT, RX, RY, RZ, CX, CZ, and CCX.

This project does not aim at implementing an especially optimized quantum statevector simulator, but rather to show the exploitation of sparse linear algebra properties to circumvent the need to construct or store any matrices larger than single qubit gates.

## Single-Qubit Gates

To optimize the simple matrix vector multiplication, we first have a look at the effect an arbitrary unitary 

```math
U=\begin{pmatrix}
   \alpha{}_{00} & \alpha{}_{01}\\
   \alpha{}_{10} & \alpha{}_{11} 
\end{pmatrix}
```

has on a given target qubit $i$ of an $n$-qubit vector:

```math
{I}^{\otimes{}2^{i-1}} \otimes{}  U \otimes{} {I}^{\otimes{}2^{n-i}}.
```

We note that ${I_k} \otimes{} {I_l} = {I_{k\cdot{}l}}$ and with that:

```math
{I_{2^{i-1}}} \otimes{}  U \otimes{} {I_{2^{n-i}}}.
```

We can then simplify to:

```math
{I_{2^{i-1}}} \otimes{}  \begin{pmatrix}
    \alpha{}_{00} \cdot{} {I_{2^{n-i}}} & \alpha{}_{01} \cdot{} {I_{2^{n-i}}} \\
    \alpha{}_{10} \cdot{} {I_{2^{n-i}}} & \alpha{}_{11} \cdot{} {I_{2^{n-i}}}
\end{pmatrix}
=
{I_{2^{i-1}}} \otimes{}  \begin{pmatrix}
    \text{diag}(\alpha{}_{00})^{2^{n-i}} & \text{diag}(\alpha{}_{01})^{2^{n-i}}  \\
    \text{diag}(\alpha{}_{10})^{2^{n-i}} & \text{diag}(\alpha{}_{11})^{2^{n-i}}
\end{pmatrix}.
```

Further put into block diagonal form:

```math
    \text{Diag}(\begin{pmatrix}
    \text{diag}(\alpha{}_{00})^{2^{n-i}} & \text{diag}(\alpha{}_{01})^{2^{n-i}}  \\
    \text{diag}(\alpha{}_{10})^{2^{n-i}} & \text{diag}(\alpha{}_{11})^{2^{n-i}}
\end{pmatrix})^{2^{i-1}}.
```

Calculations with this block diagonal matrix can then be simplified in as in the following:

1. Divide state vector into $2^{i-1}$ non-overlapping blocks of size $2^{n-i+1}$. 
2. Apply to each block the following matrix, leading to two new observations

```math
    \begin{pmatrix}
    \text{diag}(\alpha{}_{00})^{2^{n-i}} & \text{diag}(\alpha{}_{01})^{2^{n-i}}  \\
    \text{diag}(\alpha{}_{10})^{2^{n-i}} & \text{diag}(\alpha{}_{11})^{2^{n-i}}
\end{pmatrix}.
``` 


- Apply diag($\alpha{}_{00}$) and diag($\alpha{}_{10}$) to upper half of block and apply diag($\alpha{}_{01}$) and  diag($\alpha{}_{11}$) to lower half. 
- Sum appropriate applications together and write to state vector. This can be done in-place as the blocks are all non-overlapping.

## Measurement
We calculate the measurement of a given qubit by first calculating the probability of it collapsing to $\ket{0}$. We do this by applying $\ket{0}\bra{0}$ to the target qubit and then computing the norm on the resulting statevector.

For this, we can first observe that $\ket{0}\bra{0}$ only keeps those amplitudes, where the index in the statevector, when indexing in binary, has a 0 at the target qubit's place. This is because 

```math
\ket{0}\bra{0} = \begin{pmatrix} 1 & 0 \\ 0 & 0\end{pmatrix}
```

applied to $\ket{0}$ is $\ket{0}$ and to $\ket{1}$ is the 0 vector. Looking at any $I_{(2^i)} \otimes{} \ket{0}\bra{0} \otimes{} I_{(2^j)}$, we can easily see that the resulting matrix filters for all amplitudes where the target qubit is indeed inactive.

Let's say we have a 3 qubit statevector and our target qubit is qubit 2. That would mean the indices 000, 100, 001, and 101 are kept, while 010, 110, 011, and 111 are set to 0. For the first index set, we then compute the norm of the amplitudes and sum them up. 

As an algorithm, we can write:

```python
measure(statevector, target_qubit, number_qubits)
  p_0 = 0.0;
  target_mask = 1u << (number_qubits - 1 - target_qubit)
  for (i = 0; i < size(statevector); i++) 
    if ((i & target_mask) == 0) 
      p_0 += norm(statevector[i])
  p_1 = 1.0 - p_0;
```

With the probabilities of $p_0$ and $p_1$, we can compute a binomial coin flip for the measurement outcome and collapse the qubit accordingly and renormalize the statevector.

## Two Qubit Gates

With out loss of generality, a given $CX$ or $CZ$ gate can be calculated as $I \otimes{} \ket{0}\bra{0} +  U \otimes{} \ket{1} \bra{1}$, where $U$ is either a Pauli $X$ or $Z$ gate. For arbitrary distances between qubits we can write:

```math
I_{2^i} \otimes{} \ket{0} \bra{0} \otimes{} I_{2^j} \otimes{} I \otimes{} I_{2^k} + I_{2^i} \otimes{} \ket{1} \bra{1} \otimes{} I_{2^j} \otimes{} U \otimes{} I_{2^k},
```

where $i$ is the number of qubits more significant than the control qubit, $j$ the number of qubits inbetween control and target, $k$ is the number of qubits less significant than the target qubit. As an example, take a CX on $\ket{001}$, where the control qubit is the most significant and the target the least significant qubit. The formula from above then takes the explicit form of 

```math
(\ket{0}\bra{0}\otimes{}I_{2^1}\otimes{}I + \ket{1}\bra{1} \otimes{} I \otimes{} X) \ket{001} =
\begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 1 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\
\end{pmatrix}
\begin{pmatrix}
0\\
0\\
0\\
0\\
1\\
0\\
0\\
0\\
\end{pmatrix}
=
\begin{pmatrix}
0\\
0\\
0\\
0\\
0\\
1\\
0\\
0\\
\end{pmatrix}
 = \ket{101}.
```

The upper left part symbolizes the inactive case. As we have already seen in the measurement section, this means we only consider amplitudes whose binary index is 0 at the control qubit's position. We then apply the identity matrix to these amplitudes, preserving the amplitudes in the case the control qubit is 0. 

So in iterating over the statevector in the algorithm, we can already skip all those elements, since only the identity matrix is applied. This then leaves only those elements, that share an active control bit, in our example then $\ket{001}$, $\ket{101}$, $\ket{011}$, and $\ket{111}$.

Now, the one projector $\ket{1}\bra{1}$ takes all the active control qubit cases and applies the single qubit gate to the target qubit. Here, we can group those states together, that differ only in the target qubit being active or inactive. Our groups are {$\ket{001}$, $\ket{101}$} and {$\ket{011}$, $\ket{111}$}. We can now choose our representative states for each group to simply be the first element, $\ket{001}$ and $\ket{011}$. 

Keeping in line with thinking about active and inactive states, we can think about our unitary $U$ in a similar fashion, where the top row is basically the inactive case, and the bottom row the active case.

So when applying the gate matrix to our target qubit, we can take the element in the statevector at the inactive position (e.g. $\ket{011}$) and the element in the statevector at the active position (e.g. $\ket{111}$). We then add $U_{0,0}\ket{011}$ and $U_{0,1}\ket{111}$ and write the result to the element at $\ket{011}$. We then add $U_{1,0}\ket{011}$ and $U_{1,1}\ket{111}$, the active case, and write it to the element at $\ket{111}$.

This results in the following algorithm:

```python
applyControlledGate(target, control, a_00, a_01, a_10, a_11, statevector)
    control_mask = 1u << (number_qubits - 1 - control)
    target_mask = 1u << (number_qubits - 1 - target)
    for (i = 0; i < size(statevector); i++)
        if ((i & control_mask) == 0)
            continue
        if (i & target_mask)
            continue
        target_inactive_index = i
        target_active_index = i | target_mask

        inactive = statevector[target_inactive_index]
        active = statevector[target_active_index]

        statevector[target_inactive_index] = U_00 * inactive + U_01 * active
        statevector[target_active_index] = U_10 * inactive + U_11 * active
```
