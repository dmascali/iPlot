[![View iPlot on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/87709-iplot)
# iPlot - Interactive Plot
*iPlot* is a MATLAB tool for plotting 1-D or 2-D matrices interactively. Data is plotted versus row indices. 

- `iplot(X)` plots the columns in X (or their spectral amplitude) interactively, one at a time. X can be a matrix or a vector. 
- `iplot(X,Y,...)` plots the columns in X,Y..., one above the other for easy comparison.

## Functionality
*iPlot* functionalities are triggered by pressing keys (see keyboard legend below).  Main functionalities include:
- move forward or backward between the columns of the input matrix ["A" and "D" keys]
- plots the spectral amplitude of the current column [ "F" key ]
- various column ordering (*sequential*, *std +*, *std -*, *random*) ["R" key]

<img align="left"  src="https://user-images.githubusercontent.com/41842946/108685642-ab513400-74f4-11eb-8363-e2f73fd8acaa.jpg" height="300">
