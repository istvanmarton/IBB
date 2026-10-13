format long

A=load('G63.txt');
b=load('strategy_L1_G63.txt');
L1_G63=sum(abs(b' * A))

A=load('M63.txt');
b=load('strategy_L1_M63.txt');
L1_M63=sum(abs(b' * A))

A=load('M91.txt');
b=load('strategy_L1_M91.txt');
L1_M91=sum(abs(b' * A))

A=load('M97.txt');
b=load('strategy_L1_M97.txt');
L1_M97=sum(abs(b' * A))

A=load('F72.txt');
b=load('strategy_L1_F72.txt');
L1_F72=sum(abs(b' * A))

A=load('M72x360.txt');
b=load('strategy_L1_M72x360.txt');
L1_M72x360=sum(abs(b' * A))

A=load('M81x321.txt');
b=load('strategy_L1_M81x321.txt');
L1_M81x321=sum(abs(b' * A))

A=load('F90.txt');
b=load('strategy_L1_F90.txt');
L1_F90=sum(abs(b' * A))

A=load('M100.txt');
b=load('strategy_L1_M100.txt');
L1_M100=sum(abs(b' * A))

A=load('M108.txt');
b=load('strategy_L1_M108.txt');
L1_M108=sum(abs(b' * A))

A=load('M120.txt');
b=load('strategy_L1_M120.txt');
L1_M120=sum(abs(b' * A))

A=load('C120.txt');
b=load('strategy_L1_C120.txt');
L1_C120=sum(abs(b' * A))

A=load('M130.txt');
b=load('strategy_L1_M130.txt');
L1_M130=sum(abs(b' * A))
