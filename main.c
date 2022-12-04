#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tinyexpr.h"
#include "tinyexpr.c"


#define APPROXIMATE_VALUE 0.0000000000001
#define LEN_STR 10000
#define LEN_VAR_NAME 5
#define LEN_VAR_ARRAY 1000



//-------------FREE MEMORY----------------------------------
void free_array(double** array, int row){
    for(int i = 0; i < row; i++){
        free(array[i]);
    }
    free(array);
}

void free_array_char(char** array, int row){
    for(int i = 0; i < row; i++){
        free(array[i]);
    }
    free(array);
}


//------------------WORK WITH MATRIX------------------------------------
double** get_cofactor(double** matrix, int row, int col, int order){
    double** new_matrix = malloc(order*sizeof(double *));
    for(int i = 0; i < (order-1);i++){
        double* row_matrix = malloc(order*sizeof(double));
        new_matrix[i] = row_matrix;
        for(int g = 0; i < (order-1); i++){
            new_matrix[i][g] = matrix[(i + (i >= row))][(g + (g >= col))];
        }
    }
    return new_matrix;
}


double get_det(double** matrix, int order){
    double det = 0;
    if(order == 1){
        return matrix[0][0];
    } else {
        for(int i = 0; i < order; i++){
            int sign = ( i % 2 == 0) ? 1 : -1;
            det += sign*matrix[0][i]* get_det(get_cofactor(matrix, 0, i, order), (order-1));
        }
    }
    return det;
}


double** get_transposed_matrix(double** matrix, int row, int col){
    double** new_matrix = malloc(col*sizeof(double *));
    for(int i = 0; i < row; i++){
        double* row_matrix = malloc(row*sizeof(double));
        new_matrix[i] = row_matrix;
        for(int g = 0; g < col; g++){
            new_matrix[i][g] = matrix[g][i];
        }
    }
    return new_matrix;
}


double** get_inverse_matrix(double** matrix, int order){
    double** new_matrix = malloc(order*sizeof(double *));
    double det = 1.0 / (get_det(matrix, order));
    for(int i = 0; i < order; i++){
        double* row_matrix = malloc(order*sizeof(double));
        new_matrix[i] = row_matrix;
        for(int g = 0; g < order; g++){
            int sign = ( (i+g) % 2 == 0) ? 1 : -1;
            new_matrix[i][g] = sign * det * get_det(get_cofactor(matrix, i, g, order), (order-1));
        }
    }
    return get_transposed_matrix(new_matrix, order, order);
}


double** multiple_matrix(double** matrix_first, double** matrix_second, int row_first, int col_first, int row_second, int col_second){
    double** new_matrix = malloc(row_first*sizeof(double *));
    for(int i = 0; i < row_first; i++){
        double* row_matrix = malloc(col_second*sizeof(double));
        new_matrix[i] = row_matrix;
        for (int g = 0; g < col_second; g++) {
            new_matrix[i][g] = 0;
            for (int k = 0; k < col_first; k++) {
                new_matrix[i][g] += matrix_first[i][k] * matrix_second[k][g];
            }
        }
    }
    return new_matrix;
}


double** matrix_subtraction(double** matrix_1, double** matrix_2, int row, int col){
    double** new_matrix = malloc(row*sizeof(double *));
    for(int i = 0; i < row; i++) {
        double *row_matrix = malloc(col * sizeof(double));
        new_matrix[i] = row_matrix;
        for(int g = 0; g < col; g++){
            new_matrix[i][g] = matrix_1[i][g] - matrix_2[i][g];
        }
    }
    return new_matrix;
}


//------------------WORK WITH EQUATIONS---------------------------------
char** create_system(int num_equations, int num_of_vars){
    char **system_equations = malloc(num_of_vars*LEN_STR*sizeof(char));
    for(int i = 0; i < num_equations; i++){
        if(i > (num_of_vars-1)){
            char* equation = malloc(sizeof(char) * LEN_STR);
            scanf("%s", equation);
            system_equations[num_of_vars-1] = strcat(strcat(system_equations[num_of_vars-1], "+"), equation);
        } else {
            system_equations[i] = (char*)malloc(sizeof(char) * LEN_STR);
            scanf("%s", system_equations[i]);
        }
    }
    return system_equations;
}


double** form_derivative(int num_of_vars, double** array_vars, char** system_equations, te_variable variables[]){
    double** derivative_array = malloc(num_of_vars*sizeof(double*));
    int err;
    for(int i = 0; i < num_of_vars; i++){
        double* str_derivative_array = malloc(num_of_vars*sizeof(double));
        te_expr *n = te_compile(system_equations[i], variables, num_of_vars, &err);
        for(int g = 0; g < num_of_vars; g++){
            double a, b;
            a = te_eval(n);
            array_vars[g][0] = array_vars[g][0] - APPROXIMATE_VALUE;
            b = te_eval(n);
            array_vars[g][0] = array_vars[g][0] + APPROXIMATE_VALUE;
            str_derivative_array[g] = (a-b)/(APPROXIMATE_VALUE);
        }
        te_free(n);
        derivative_array[i] = str_derivative_array;
    }
    return derivative_array;
}


double** form_func_arr(int num_of_vars, double** array_vars, char** system_equations, te_variable variables[]){
    double** func_array = malloc(num_of_vars*sizeof(double*));
    int err;
    for(int i = 0; i < num_of_vars; i++){
        double* str_derivative_array = malloc(1*sizeof(double));
        te_expr *n = te_compile(system_equations[i], variables, num_of_vars, &err);
        str_derivative_array[0] = te_eval(n);
        func_array[i] = str_derivative_array;
        te_free(n);
    }
    return func_array;
}


int main() {
    int num_equations = 0, num_of_vars = 0, err = 0;
    double guess = 3;


    printf("Please, enter number of variables from the system of equations:\n");
    scanf("%d", &num_of_vars);


    double** array_vars = malloc(num_of_vars*sizeof(double*));
    te_variable variables[LEN_VAR_ARRAY] = {};
    char** name_of_vars = malloc(num_of_vars*sizeof(char *));


    printf("Please, enter variables from the system of equations separated by a space:\nExample: x y w z\n");
    for(int i = 0; i < num_of_vars; i++){
        double* array_len_vars = malloc(1*sizeof(double));
        array_vars[i] = array_len_vars;
        char* name_var = malloc(LEN_VAR_NAME*sizeof(char));
        scanf("%s", name_var);
        name_of_vars[i] = name_var;
        array_vars[i][0] = guess;
        te_variable variable = {name_of_vars[i], &array_vars[i][0]};
        variables[i] = variable;
    }


    printf("Enter the hypothetical root of the equations:\n");
    scanf("%lf", &guess);
    printf("Enter number of equations:\n");
    scanf("%d", &num_equations);
    if(num_of_vars < num_equations){
        printf("WARNING!!!\nIf there are more equations than the number of variables, then the program may give an incorrect result!\n");
    } else if(num_of_vars > num_equations) {
        printf("WARNING!!!\nNumber of equations < number of variables\n");
        return 0;
    }
    printf("Enter equations:\n");
    printf("Examples:\nx^2+y^2-2\nx-y+3\nsqrt(x+y)+2\n");
    printf("Enter equations:\n");
    char** system_equations = create_system(num_equations, num_of_vars);


    for(int i = 0; i < 1000; i++){
        double** derivative_array = get_inverse_matrix(form_derivative(num_of_vars, array_vars, system_equations, variables), num_of_vars);
        double** func_array = form_func_arr(num_of_vars, array_vars, system_equations, variables);
        double** mult_array = multiple_matrix(derivative_array, func_array, num_of_vars, num_of_vars, num_of_vars, 1);
        array_vars = matrix_subtraction(array_vars, mult_array, num_of_vars, 1);
        for(int i = 0; i < num_of_vars; i++){
            te_variable variable = {name_of_vars[i], &array_vars[i][0]};
            variables[i] = variable;
        }
        free_array(derivative_array, num_of_vars);
        free_array(func_array, num_of_vars);
        free_array(mult_array, num_of_vars);
    }


    printf("\n");
    for(int i = 0; i < num_of_vars; i++){
        printf("%s = %lf\n", name_of_vars[i], array_vars[i][0]);
    }


    free_array(array_vars, num_of_vars);
    free_array_char(system_equations, num_of_vars);
    free_array_char(name_of_vars, num_of_vars);
    return 0;
}