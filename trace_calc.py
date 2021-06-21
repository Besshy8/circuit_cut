import numpy as np
## from qiskit_textbook.tools import array_to_latex

class Tr:
    """
    This is a Python Class to compute traces and partial traces of (4 times 4) matrices.

    # Example

    If you want to calculate the matrix M using this class, you can call the method in the class as follows

    ```python
        ans = Tr(M).p("b")
    ```

    The argument "p" can specify how to take the partial trace. If it is "a", the partial trace is calculated for the first qubit,
    and if it is "b", the partial trace is calculated for the second qubit.
    If no argument is specified, a simple trace will be taken.

    if you have qiskit_textbook.tools , you can display partial trace using tex format.

    # Arguments
        M: (4 times 4) matrices which you want to calc traces or partial traces

    """

    def __init__(self, M):
        self.M = M
        self.div_mat = self._division_mat(self.M)

    def p(self, flag=None):
        if flag == None:
            return np.trace(self.M)
        else:
            self.flag = flag
            self.partial_tr = np.zeros((2, 2), dtype=np.complex)
            c = 0
            for k in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                a = self.div_mat[c]
                c += 1
                for i in [0, 1]:
                    for j in [0, 1]:
                        if self.flag == "a": self.partial_tr += a[i, j] * np.trace(
                            np.dot(self._ket(k[0]), self._bra(k[1]))) * np.dot(self._ket(i), self._bra(j))
                        if self.flag == "b": self.partial_tr += a[i, j] * np.dot(self._ket(k[0]),
                                                                                 self._bra(k[1])) * np.trace(
                            np.dot(self._ket(i), self._bra(j)))
        ## return array_to_latex(self.partial_tr)  ## If you want to use the tex format, uncomment it!
        return self.partial_tr

    def _bra(self, i):
        if i == 0:
            ans = np.array([[1, 0]])
        elif i == 1:
            ans = np.array([[0, 1]])
        return ans

    def _ket(self, j):
        if j == 0:
            ans = np.array([[1], [0]])
        elif j == 1:
            ans = np.array([[0], [1]])
        return ans

    def _division_mat(self, M):
        self.upper_left = M[0:2, 0:2]
        self.lower_right = M[0:2, 2:4]
        self.upper_right = M[2:4, 0:2]
        self.lower_left = M[2:4, 2:4]
        return (self.upper_left, self.lower_right, self.upper_right, self.lower_left)


if __name__ == '__main__':
    M = np.array(np.arange(16).reshape((4, 4)))
    print("Matrix: ")
    print(M)
    print("Trace: ", Tr(M).p())
    print("Trace_A: ")
    print(Tr(M).p("a"))
    print("Trace_B: ")
    print(Tr(M).p("b"))