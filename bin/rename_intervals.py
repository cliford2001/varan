import glob
import os

# Obtener la lista de archivos de intervalos en el directorio de salida
intervals = sorted(glob.glob("*_split/*/*.interval_list"))

# Iterar sobre cada archivo de intervalos y renombrarlo
for i, interval in enumerate(intervals):
    # Obtener la ruta y el nombre del archivo
    directory, filename = os.path.split(interval)
    
    # Construir el nuevo nombre de archivo con un prefijo numérico para evitar conflictos de nombres
    new_filename = os.path.join(directory, f"{i + 1}_{filename}")
    
    try:
        # Renombrar el archivo
        os.rename(interval, new_filename)
        print(f"Renombrado '{filename}' a '{os.path.basename(new_filename)}'")
    except Exception as e:
        print(f"Error al renombrar '{filename}': {e}")

