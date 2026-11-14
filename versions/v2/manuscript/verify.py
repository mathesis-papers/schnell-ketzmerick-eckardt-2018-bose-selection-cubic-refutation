import itertools
import json
from fractions import Fraction
from math import comb


def determinant(matrix):
    if not matrix:
        return Fraction(1)
    result = Fraction(0)
    for column, entry in enumerate(matrix[0]):
        minor = [row[:column] + row[column + 1:] for row in matrix[1:]]
        result += (-1) ** column * entry * determinant(minor)
    return result


def pfaffian_vector(matrix):
    result = []
    for removed in range(5):
        first, second, third, fourth = [index for index in range(5) if index != removed]
        result.append((-1) ** removed * (
            matrix[first][second] * matrix[third][fourth]
            - matrix[first][third] * matrix[second][fourth]
            + matrix[first][fourth] * matrix[second][third]))
    return result


def add_polynomials(*polynomials):
    result = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            result[monomial] = result.get(monomial, 0) + coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def multiply_polynomials(first, second):
    result = {}
    for first_term, first_coefficient in first.items():
        for second_term, second_coefficient in second.items():
            monomial = tuple(sorted(first_term + second_term))
            result[monomial] = result.get(monomial, 0) + first_coefficient * second_coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def symbolic_kernel_check():
    matrix = [[{} for column in range(5)] for row in range(5)]
    for row, column in itertools.combinations(range(5), 2):
        matrix[row][column] = {((row, column),): 1}
        matrix[column][row] = {((row, column),): -1}
    pfaffians = []
    for removed in range(5):
        first, second, third, fourth = [index for index in range(5) if index != removed]
        terms = [multiply_polynomials(matrix[first][second], matrix[third][fourth]),
                 multiply_polynomials(matrix[first][third], matrix[second][fourth]),
                 multiply_polynomials(matrix[first][fourth], matrix[second][third])]
        terms = [{monomial: coefficient * (-1) ** (removed + position % 2)
                  for monomial, coefficient in term.items()}
                 for position, term in enumerate(terms)]
        pfaffians.append(add_polynomials(*terms))
    for row in range(5):
        polynomial = add_polynomials(*[
            multiply_polynomials(matrix[row][column], pfaffians[column])
            for column in range(5)])
        assert not polynomial, (row, polynomial)


def higher_level_checks(first, second):
    records = []
    for dimension in (6, 7, 8):
        added = dimension - 5
        strength = Fraction(1, 1000)
        first_full = [(value / (1 + 2 * added * strength))
                      for value in first + [2 * strength] * added]
        second_full = [(value / (1 + added * strength))
                       for value in second + [strength] * added]
        assert sum(first_full) == sum(second_full) == 1
        matrix = [[(column - row) ** 3
                   * (first_full[row] * first_full[column]
                      - second_full[row] * second_full[column])
                   for column in range(dimension)] for row in range(dimension)]
        principal = [row[:5] for row in matrix[:5]]
        pfaffians = pfaffian_vector(principal)
        assert all(value > 0 for value in pfaffians)
        occupation = [value / sum(pfaffians) for value in pfaffians] + [Fraction(0)] * added
        drift = [sum(entry * weight for entry, weight in zip(row, occupation)) for row in matrix]
        assert drift[:5] == [0] * 5
        assert all(value < 0 for value in drift[5:])
        assert determinant([row[1:] for row in principal[1:]]) > 0
        records.append({'dimension': dimension, 'added_weight_parameter': str(strength),
                        'selected_modes': [1, 2, 3, 4, 5],
                        'outside_drifts': [str(value) for value in drift[5:]]})
    return records


def monomial_factorization_checks(first, second):
    records = []
    energies = [1, 2, 3, 4, 5]
    for exponent in (1, 3, 5):
        feature_count = 2 * (exponent + 1)
        features = [[weight[row] * energies[row] ** degree
                     for weight in (first, second) for degree in range(exponent + 1)]
                    for row in range(5)]
        form = [[0] * feature_count for row in range(feature_count)]
        for bath, sign in enumerate((1, -1)):
            for degree in range(exponent + 1):
                form[bath * (exponent + 1) + degree][bath * (exponent + 1) + exponent - degree] = (
                    sign * (-1) ** degree * comb(exponent, degree))
        assert all(form[row][column] == -form[column][row]
                   for row in range(feature_count) for column in range(feature_count))
        for row in range(5):
            for column in range(5):
                reconstructed = sum(features[row][left] * form[left][right] * features[column][right]
                                    for left in range(feature_count) for right in range(feature_count))
                direct = (energies[column] - energies[row]) ** exponent * (
                    first[row] * first[column] - second[row] * second[column])
                assert reconstructed == direct
        records.append({'exponent': exponent, 'bath_count': 2, 'feature_columns': feature_count,
                        'entry_identities': 25, 'factorization': 'A = V H V^T'})
    return records


def main():
    energies = [1, 2, 3, 4, 5]
    first = [Fraction(value, 21) for value in [5, 4, 4, 4, 4]]
    second = [Fraction(value, 21) for value in [6, 6, 3, 3, 3]]
    assert sum(first) == sum(second) == 1
    assert all(value > 0 for value in first + second)
    matrix = [[(energies[column] - energies[row]) ** 3
               * (first[row] * first[column] - second[row] * second[column])
               for column in range(5)] for row in range(5)]
    integer_matrix = [[0, -16, 16, 54, 128], [16, 0, -2, -16, -54],
                      [-16, 2, 0, 7, 56], [-54, 16, -7, 0, 7],
                      [-128, 54, -56, -7, 0]]
    assert matrix == [[Fraction(value, 441) for value in row] for row in integer_matrix]
    assert all(matrix[row][column] == -matrix[column][row]
               for row in range(5) for column in range(5))
    occupation = [Fraction(value, 100) for value in [14, 56, 21, 8, 1]]
    assert sum(occupation) == 1 and all(value > 0 for value in occupation)
    residual = [sum(entry * weight for entry, weight in zip(row, occupation)) for row in matrix]
    assert residual == [0] * 5
    pfaffians = pfaffian_vector(integer_matrix)
    assert pfaffians == [504, 2016, 756, 288, 36]
    assert pfaffian_vector(matrix) == [Fraction(value, 441 ** 2) for value in pfaffians]
    assert [Fraction(value, sum(pfaffians)) for value in pfaffians] == occupation
    principal = determinant([row[1:] for row in integer_matrix[1:]])
    assert principal == 254016 == pfaffians[0] ** 2
    assert determinant(integer_matrix) == 0
    symbolic_kernel_check()
    higher_levels = higher_level_checks(first, second)
    factorizations = monomial_factorization_checks(first, second)
    print(json.dumps({
        'arithmetic': 'exact rational and integer polynomial coefficients',
        'coupling_weights': {'f': [str(value) for value in first],
                             'g': [str(value) for value in second]},
        'integer_matrix': integer_matrix,
        'matrix_scale': '1/441',
        'occupation_fractions': [str(value) for value in occupation],
        'kernel_residual': [str(value) for value in residual],
        'signed_pfaffians_integer_matrix': pfaffians,
        'principal_minor_determinant': int(principal),
        'matrix_rank': 4,
        'polynomial_kernel_identities': 5,
        'higher_level_checks': higher_levels,
        'monomial_factorizations': factorizations,
        'checks': 'pass'
    }, indent=2))


if __name__ == '__main__':
    main()
