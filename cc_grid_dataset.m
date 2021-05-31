function [x,w]=cc_grid_dataset ( m, n_1d )

%*****************************************************************************80
%
%% TENSOR_GRID_DISPLAY displays a tensor product grid.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    30 December 2011
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, the spatial dimension.
%    Only M = 1, 2, and 3 are handled by this program.
%
%    Input, integer N_1D(M), the order of the 1D rules.
%    If N_1D is a scalar, it will be extended to a vector of length M,
%    with all entries equal.
%
%
%  Local Parameters:
%
%    Local, integer N, the number of points in the tensor product grid.
%
%    Local, real X(M,N), the points in the tensor product grid.
%
 
%  Make sure we have the arguments.
%
  if ( nargin < 1 )
    m = input ( 'Enter the spatial dimension M, 1, 2 or 3:  ' );
  elseif ( ischar ( m ) )
    m = str2num ( m );
  end

  if ( nargin < 2 )
    n_1d = input ( 'Enter the number of points in the 1D grid, N1D:  ' );
  elseif ( ischar ( n_1d ) )
    n_1d = str2num ( n_1d );
  end

  if ( 1 < m && length ( n_1d ) == 1 )
    n_1d = n_1d * ones ( m, 1 );
  end

  [x, w] = tensor_product ( m, n_1d);
  point_num=length(w);
  level_max=max(n_1d);

%
%  Write the rule to files.
%
  w_filename = sprintf ( 'ccfull_d%d_level%d_w.txt', m, level_max );
  x_filename = sprintf ( 'ccfull_d%d_level%d_x.txt', m, level_max );

  %fprintf ( 1, '  Creating W file = "%s".\n', w_filename );

  r8mat_write ( w_filename, 1, point_num, w );

  %fprintf ( 1, '  Creating X file = "%s".\n', x_filename );

  r8mat_write ( x_filename, m, point_num, x );
  


  return
end

function [ x, w ] = cc_compute ( n )

%*****************************************************************************80
%
%% CC_COMPUTE computes a Clenshaw Curtis quadrature rule.
%
%  Discussion:
%
%    Our convention is that the abscissas are numbered from left to right.
%
%    The rule is defined on [-1,1].
%
%    The integral to approximate:
%
%      Integral ( -1 <= X <= 1 ) F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    15 February 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the order of the rule.
%    1 <= N.
%
%    Output, real X(N), the abscissas.
%
%    Output, real W(N), the weights.
%
  if ( n < 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'CC_COMPUTE - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of N = %d\n', n );
    error ( 'CC_COMPUTE - Fatal error!' );
  end

  w = zeros ( n, 1 );
  x = zeros ( n, 1 );

  if ( n == 1 )
    x(1) = 0.0;
    w(1) = 2.0;
    return
  end

  for i = 1 : n
    x(i) = cos ( ( n - i ) * pi / ( n - 1 ) );
  end

  x(1) = -1.0;
  if ( mod ( n, 2 ) == 1 )
    x((n+1)/2) = 0.0;
  end
  x(n) = +1.0;

  for i = 1 : n

    theta = ( i - 1 ) * pi / ( n - 1 );

    w(i) = 1.0;

    for j = 1 : ( n - 1 ) / 2

      if ( 2 * j == ( n - 1 ) )
        b = 1.0;
      else
        b = 2.0;
      end

      w(i) = w(i) - b * cos ( 2.0 * j * theta ) / ( 4 * j * j - 1 );

    end

  end

  w(1)     =       w(1)     / ( n - 1 );
  w(2:n-1) = 2.0 * w(2:n-1) / ( n - 1 );
  w(n)     =       w(n)     / ( n - 1 );

  return
end

function [grid_x, grid_w] = tensor_product ( dim_num, order )

%*****************************************************************************80
%
%% TENSOR_PRODUCT computes the points of a tensor product grid.
%
%  Discussion:
%
%    If the coordinates of each factor rule are given in lexicographic order,
%    then the DIM_NUM-dimensional points produced by this function will
%    also be in lexicographic order.
%
%  Example:
%
%    DIM_NUM = 3
%    ORDER = [ 2, 3, 2 ]
%    RULE = { [1], [2,3], [4,5,6], ... } (just for clarity)
%
%    GRID_NUM = 12
%    GRID_X =
%      2  2  2  2  2  2  3  3  3  3  3  3
%      4  4  5  5  6  6  4  4  5  5  6  6
%      2  3  2  3  2  3  2  3  2  3  2  3
%
%    Read the columns of GRID_X to see the 12 vectors created.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    09 October 2011
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer DIM_NUM, the spatial dimension.
%
%    Input, integer ORDER(DIM_NUM), the order of the rule in each dimension.
%
%    Input, integer RULE, the index of the rule used in all dimensions.
%
%    Output, real GRID_X(DIM_NUM,GRID_NUM), the points of the grid.
%
  grid_num = prod ( order(1:dim_num) );

  grid_x = zeros ( dim_num, grid_num );
  grid_w = zeros ( 1, grid_num );

  chunk_num = 1;
  rep_num = grid_num;

  for dim = 1 : dim_num

    [x,w] = cc_compute ( order(dim) );
    x = x ( : );
    w = w ( : );
    
    chunk_num = chunk_num * order(dim);
    rep_num = rep_num / order(dim);

    i = 1;
    j = 1;
    for chunks = 1 : chunk_num
      for reps = 1 : rep_num
        grid_x(dim,j) = x(i);
        if dim==1 
            grid_w(j) = w(i);
        end
        j = j + 1;
      end
      i = i + 1;
      if ( order(dim) < i )
        i = 1;
      end
    end

  end
  grid_w=grid_w./sum(grid_w); %normalize

  return
end

function r8mat_write ( output_filename, m, n, table )

%*****************************************************************************80
%
%% R8MAT_WRITE writes an R8MAT file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 August 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string OUTPUT_FILENAME, the output filename.
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of points.
%
%    Input, real TABLE(M,N), the points.
%

%
%  Open the file.
%
  output_unit = fopen ( output_filename, 'wt' );

  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8MAT_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'R8MAT_WRITE - Error!' );
  end
%
%  Write the data.
%
%  For smaller data files, and less precision, try:
%
%     fprintf ( output_unit, '  %14.6f', table(i,j) );
%
  for j = 1 : n
    for i = 1 : m
      fprintf ( output_unit, '  %24.16f', table(i,j) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  Close the file.
%
  fclose ( output_unit );

  return
end