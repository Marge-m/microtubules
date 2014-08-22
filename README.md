microtubules
============

Кинетическая модель сборки-разборки микротрубочек
MTevolve.m - запускает event.m (преобразования матриц с димерами и с латеральными связями)

MTevolve(N, M), где N - максимальная длина микротрубочки (задавать данный параметр следует исходя из того, что за 1 секунду прирастает до 8 димеров), M - количество операций (где операция одно из событий присоединение димера/отрыв димера/образование латеральной связи/разрыв латеральной связи/гидролиз GTP, присоединенной к тубулину). 1000 операций соответствуют ~0.75 сек реального времени.

Модель составлена по статьям

1. The mechanisms of microtubule catastrophe and rescue: implications from analysis of a dimer-scale computational model
Gennady Margolin, Ivan V. Gregoretti, Trevor M. Cickovski, Chunlei Li, Wei Shi, Mark S. Alber and Holly V. Goodsonb, 2012.
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3279392/

2. Estimates of lateral and longitudinal bond energies within the microtubule lattice
Vincent VanBuren, David J. Odde, and Lynne Cassimeris, 2002
http://www.pnas.org/content/99/9/6035.full



Развернутая микротрубочка в фазу роста:

![Развернутая микротрубочка в фазу роста](http://i.imgur.com/6YOjuYR.gif)
