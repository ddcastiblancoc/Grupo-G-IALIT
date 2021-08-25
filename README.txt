============================================================
README:


Para realizar la simulación de los objetivos planteados se tienen los siguientes códigos:

1) maxwell_juttner_algorithms.py
2) maxwell_juttner_results.py
3) PIC_algorithms.py
4) PIC_initialization.py
5) PIC_loop.py
6) PIC_results.py

============================================================

Para correrlos, todos los archivos deben estar ubicados en una misma carpeta.

1) y 3) sólo es necesarios correrlos una sola vez para cargar las funciones a utilizar.

Sus funciones son importadas a los otros códigos con import maxwell_juttner_algorithms as mj y import PIC_algorithms as PIC.

============================================================

RESPECTO A LOS ARCHIVOS MAXWELL_JUTTNER:

2) permite recrear las gráficas de las distribuciones estacionarias, desplazadas y la comparación entre los tres algoritmos que fueron mostradas en el informe.

Puede variar el número total de números muestreados (no recomendable hacerlo más de 10 millones). Para la gráfica de las distribuciones por Sobol puede variar las temperaturas theta (mantenga la misma cantidad de valores) y velocidades de boost beta. Puede cambiar entre los tres algoritmos (rechazo, Sobol, transformada inversa) en la línea 57 según el nombre de las funciones en 1).

Debe advertirse que para valores de theta < 1, Sobol deja de aceptar valores y no funciona adecuadamente, por lo que NO se recomienda probar con esos valores.

Para la gráfica de error puede variar theta y beta a gusto.

============================================================

RESPECTO A LOS ARCHIVOS PIC:

>Si se quiere compilar de cero:

Estos deben ejecutarse en orden: 4) -> 5) -> 6).



4) calcula posiciones, velocidades de partículas del medio y campos eléctrico y magnético iniciales del sistema.

Compilar 4) exportará en la carpeta los archivos:
'particles_i_IGM.csv'
'particles_e_IGM.csv'
'E_x.csv'
'E_y.csv'
'B_z.csv'

Esto se hace para que los valores de esta inicialización quede guardada de forma segura y puedan emplearse para la segunda etapa del algoritmo PIC.

Puede cambiar en este el tamaño de grilla y el número de partículas del medio (preferiblemente iguales entre iones y electrones).



5) implementa las 5 etapas del ciclo PIC y da los resultados de las posiciones, velocidades de partículas del medio/jet y campos eléctrico y magnético finales del sistema.

Compilar 5) tomará los archivos .csv y exportará en la carpeta los archivos:
'particles_i_IGM_final.csv'
'particles_e_IGM_final.csv'
'particles_i_jet_final.csv'
'particles_e_jet_final.csv'
'E_x_final.csv'
'E_y_final.csv'
'B_z_final.csv'

Si se cambió el tamaño de la grilla en la inicialización, debe cambiarse acorde también acá. Puede cambiarse a gusto el ancho del jet, el factor de Lorentz de inyección, el número de partículas introducidas por el jet en cada instante de tiempo o el número total de pasos de tiempo.

Debe tenerse en cuenta que dada la complejidad de operaciones a realizar el código puede tardarse más tiempo entre mayor sea el tamaño de la grilla, el número de partículas en el medio o el tiempo total de pasos temporales. Para los valores que vienen por defecto, el código puede demorarse de 7-8 horas o más, dependiendo del computador. Si se quiere disminuir el tiempo, disminuya estos valores.



6) calcula gráficas basadas en las posiciones en el instante final y el campo magnético final.

Compilar 6) tomará los archivos .csv y hará las respectivas gráficas con los archivos corriendo todo de cero.



>Si no se quiere compilar de cero (OPCIONAL):

En las carpetas "dm" y "fm" vienen los datos finales ya compilados por nosotros después de varias horas de un caso débilmente magnetizado y otro fuertemente magnetizado, que puede usar directamente con 6) para obtener las gráficas presentadas en el informe.

Recuerde comentar las secciones que no va a graficar.