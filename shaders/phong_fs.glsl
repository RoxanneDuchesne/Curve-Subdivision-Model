#version 330 core
out vec4 fragColor;

in Data {
        vec3 position;
        vec3 normal;
} data;

uniform vec3 lightPosition;
uniform vec3 viewPosition;
uniform vec3 lightColor = vec3(0.3, 0, 0.9);

void main()
{
    float ambient_str = 0.2;
    float diffuse_str = 0.5;
    float specular_str = 0.5;

    //ambient lighting calcs
    vec3 ambient_vec = ambient_str * lightColor;

    //diffuse lighting calcs
    vec3 light_dir = normalize(data.position - lightPosition);
    float diffuse = max(dot(normalize(data.normal), light_dir), 0.0);
    vec3 diffuse_vec = (diffuse_str * diffuse) * lightColor;

    //specular lighting calcs
    vec3 view_dir = normalize(data.position - (-viewPosition));
    vec3 reflect_dir = reflect(light_dir, normalize(data.normal));
    float specular = pow(max(dot(reflect_dir, view_dir), 0.0), 32);
    vec3 specular_vec = (specular_str * specular) * lightColor;

    //color is combo of ambient and diffuse lighting
    fragColor = vec4(ambient_vec + diffuse_vec + specular_vec, 1.0);
   // fragColor = vec4(ambient_vec, 1.0);
    //fragColor = vec4(data.normal, 1.0);
}
