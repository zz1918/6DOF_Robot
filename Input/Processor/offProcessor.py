def convert_txt_to_off(input_filename, output_filename):
    with open(input_filename, 'r') as infile:
        lines = infile.readlines()
        num_vertices = int(lines[0].strip())
        vertices = lines[1:num_vertices + 1]
        num_faces = int(lines[num_vertices + 1].strip())
        faces = lines[num_vertices + 2:num_vertices + 2 + num_faces]

    with open(output_filename, 'w') as outfile:
        outfile.write("OFF\n")
        outfile.write(f"{num_vertices} {num_faces}\n")
        
        # Write vertices
        for vertex in vertices:
            vertex = vertex.strip()
            if vertex:
                outfile.write(vertex + "\n")
        
        # Write faces
        for face in faces:
            indices = face.strip().split()
            if len(indices) > 1:
                outfile.write(f"{len(indices) - 1} " + " ".join(indices[1:]) + "\n")

input_filename = 'output-tmp-py.txt'
output_filename = 'output.off'

convert_txt_to_off(input_filename, output_filename)
print(f"Conversion complete. The OFF file is saved as {output_filename}.")
