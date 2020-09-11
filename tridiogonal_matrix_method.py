/*
Данная программа реализует численный метод подгонки для
трёхдиагональных матриц с диагональным преобладанием
*/

class LinearEquationSystem:

    def __init__(self, size, Matrix, f):
        
        self.size = size
        self._A = Matrix
        self._f = f
        self._x = []
        self.alphas = [0]
        self.betas = [0]

    def find_coefs(self):

        self.alphas.append(-self._A[0][1] / self._A[0][0])
        self.betas.append(self._f[0] / self._A[0][0])

        i = 1
        while (i != (self.size - 1)):

            tmp = self._A[i][i - 1] * self.alphas[i] + self._A[i][i]

            alpha = -self._A[i][i + 1] / tmp
            beta = (self._f[i] - self._A[i][i - 1] * self.betas[i]) / tmp
            self.alphas.append(alpha)
            self.betas.append(beta)
            i += 1 

    def solve(self):
        
        self.find_coefs()
        tmp = self._A[self.size - 1][self.size - 2] * self.alphas[self.size - 1] + self._A[self.size - 1][self.size - 1]
        x = (self._f[self.size - 1] - self._A[self.size - 1][self.size - 2] \
             * self.betas[self.size - 1]) / tmp
        self._x.append(x)

        N = self.size - 1
        i = 0
        while (N != 0):

            self._x.append(self.alphas[N] * self._x[i] + self.betas[N])
            N -= 1
            i += 1

def main():
   
    print("Input size of matrix:")
    n = int(input())
    print("Input tridiagonal matrix:")

    matrix = []

    for _ in range(n):
        
        tmp = list(map(int, input().split()))
        matrix.append(tmp)
    
    print("Input free vector:")

    f = list(map(int, input().split()))

    flag = matrix[0][0] > matrix[0][1]
    for i in range(1, n):
        assert flag,  "Error. Your matrix of coefficients without diagonal predominance."
        if i == n - 1:
            flag = matrix[i][i] > matrix[i][i - 1]
        else:
            flag = matrix[i][i] > matrix[i][i - 1] + matrix[i][i + 1]

    sys = LinearEquationSystem(n, matrix, f)

    sys.solve()
    print("The solution of the system:")
    print(*sys._x)

main()
