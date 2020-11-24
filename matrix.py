from typing import NewType
from math import sin, cos, radians, sqrt
import pprint
Matrix = NewType("Matrix", list)

class Matrix:
    """An implementation of a Matrix object, with support for several basic operations"""
    def __init__(self, matrix=None):
        if matrix == None:
            self.matrix = []
            self.identity = []
        else:
            cols = len(matrix[0])
            rows = 0
            for index_row, row in enumerate(matrix):
                rows += 1
                assert cols == len(row), "Matrix has rows of different sizes!"
            self.dimensions = (cols, rows)
            self.matrix = matrix
            self.identity = [[1 if index_col == index_row
                              else 0 for index_col, col in enumerate(row)]
                             for index_row, row in enumerate(matrix)]

    def identity():
        return Matrix(self.identity)

    def _eq_sizes(self, other_matrix):
        if self.dimensions != other_matrix.dimensions:
            return False
        else:
            return True

    def __setitem__(self, key, value):
        assert len(self.matrix[key]) == len(value), "Length of input exceeds matrix dimensions (__setitem__)"
        self.matrix[key] = value

    def __add__(self, maddened): #maddened, matrix addend
        if isinstance(maddened, Matrix):
            assert self._eq_sizes(maddened), "Matrix must have the same dimensions."
            return Matrix([[value + maddened[index_row][index_col]
                        for index_col, value in enumerate(row)]
                       for index_row, row in enumerate(self.matrix)])
        if (isinstance(maddened, int)
            or isinstance(maddened, float)
            or isinstance(maddened, complex)):
            identadd = maddened * Matrix(self.identity)
            return [[self.matrix[index_row][index_col]
                     + identadd[index_row][index_col]
                     for index_col, col in enumerate(row)]
                    for index_row, row in enumerate(self.matrix)]

    def __mul__(self, multiplicand):
        #TODO: make multiply function
        if (isinstance(multiplicand, int)
            or isinstance(multiplicand, float)
            or isinstance(multiplicand, complex)):
            return Matrix([[value * multiplicand for value in row] for row in self.matrix])

        if (isinstance(multiplicand, Matrix)):
            assert ((cols := self.dimensions[0]) ==
                    (rows := multiplicand.dimensions[1])), "Invalid matrix dimension"
            return Matrix([[sum(value * multiplicand[index_col][i]
                                 for index_col,value in enumerate(row))
                            for i in range(cols)]
                           for row in self.matrix])
        if (isinstance(multiplicand, Vector)):
            assert((cols := self.dimensions[0]) == (rows := multiplicand.length)), "For a dot product, the vector must have the same number of elements as there are columns in the matrix"
            return Vector([sum([value * multiplicand[index] for index, value in enumerate(row)]) for row in self.matrix])
        else:
            return NotImplemented
    
    def __rmul__(self, multiplicand):
        #TODO: make multiply function
        if (isinstance(multiplicand, int)
            or isinstance(multiplicand, float)
            or isinstance(multiplicand, complex)):
            return Matrix([[value * multiplicand for value in row] for row in self.matrix])
        if (isinstance(multiplicand, Matrix)):
            assert ((cols := multiplicand.dimensions[0])
                    == (rows := self.dimensions[1])), "Invalid matrix dimensions"

            return Matrix([[sum(value * self.matrix[index_col][i]
                                  for index_col,value in enumerate(row))
                             for i in range(cols)]
                            for row in multiplicand.matrix])

    def transpose(self):
        if self.dimensions[0] <= 1 or self.dimensions[1] <=1:
            return self
        return Matrix([[self.matrix[i][col_index] for i in range(self.dimensions[1]) ]
                       for col_index, col in enumerate(self.matrix[0])])
    t = transpose

    def determinant(self):
        det = 0
        width, height = self.dimensions
        assert width == height, "A non-square matrix cannot have a determinant."
        if width == 1:
            return self.matrix[0][0]
        if width == 2:
            return (self.matrix[0][0] * self.matrix[1][1]
                    - self.matrix[0][1] * self.matrix[1][0])
        if width >= 3:
            row = self.matrix[0]
            for index, coefficient in enumerate(row):
                if index % 2 == 1:
                    submatrix = [[value for i, value in enumerate(row)
                                  if i != index] for row in self.matrix[1:]]
                    det -= coefficient * Matrix(submatrix).determinant()
                else:
                    submatrix = [[value for i, value in enumerate(row)
                                  if i != index] for row in self.matrix[1:]]
                    det += coefficient * Matrix(submatrix).determinant()
        return det
    det = determinant

    def inverse(self):
        width, height = self.dimensions
        cofactor_matrix = []
        for row_num, row in enumerate(self.matrix):
            cur_row = []
            for col_num, cell in enumerate(row):
                minor = [[cl for j_cl, cl in enumerate(rw) if j_cl != col_num]
                         for i_rw, rw in enumerate(self.matrix) if i_rw != row_num]
                cell = (Matrix(minor).determinant())
                if (row_num + col_num) % 2 == 0:
                    cur_row.append(cell)
                else:
                    cur_row.append(-cell)
            cofactor_matrix.append(cur_row)
        adjugate = Matrix(cofactor_matrix).transpose()
        return 1/self.determinant() * adjugate
    def __getitem__(self, row_key):
        return Vector(self.matrix[row_key])

    def __eq__(self, matrixcomp: Matrix) -> bool:
        if isinstance(matrixcomp, Matrix):
            if not self._eq_sizes(matrixcomp):
                return False
            for index_row, row in enumerate(self.matrix):
                for index_col, value in enumerate(row):
                    if value != matrixcomp[index_row][index_col]:
                        return False
        else:
            return False

        return True

    def __str__(self):
        return str(self.matrix)
    __repr__ = __str__

class Vector:
    """An implementation of a vector with some elementary operations."""
    def __init__(self,array):
        self.array = array
        self.length = len(array)

    def __str__(self):
        return str(self.array)
    __repr__ = __str__

    def __add__(self,addened):
        if isinstance(addened, Vector):
            assert addened.length == self.length, "Length of rows must be the same in an elementary row operation"
            add = lambda x, y: x + y
            return list(map(add,self.array, addened))
        if (isinstance(addened, int)
            or isinstance(addened, float)
            or isinstance(addened, complex)):
            return [x + addened for x in self.array]

    def __getitem__(self, index):
        return self.array[index]

    def __setitem__(self,index,value):
        self.array[index] = value

    def dot_product(self, multiplicand):
        assert(self.length == multiplicand.length), "Dot product requires that both vectors be the same length."
        return sum(val * multiplicand.array[index] for index, val in enumerate(self.array))
    def magnitude(self):
        sum = 0
        for i in self.array:
            sum += i
        return (sqrt(sum))

    def transpose(self):
        return Matrix([self.array]).transpose()
    def cross_multiply(self, multiplicand):
        assert multiplicand.length == self.length and multiplicand.length == 3, "Length of vectors must be the same when cross multiplying"
        return Vector([self.array[1] * multiplicand.array[2] - self.array[2] * multiplicand.array[1],
                self.array[2] * multiplicand.array[0] - self.array[0] * multiplicand.array[2],
                self.array[0] * multiplicand.array[1] - self.array[1] * multiplicand.array[0]])


def Rotation_Matrix(alpha):
    """Returns the matrix [[cos(alpha), sin(alpha), 0],[sin(alpha), cos(alpha), 0],[0, 0, 1]]. Alpha should be in DEGREES, not radians"""
    return Matrix([[cos(radians(alpha)), -sin(radians(alpha)), 0],
                   [sin(radians(alpha)), cos(radians(alpha)),0],
                   [0, 0, 1]])
