#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>
#include "graphics/vec.h"
#include "graphics/ivec.h"
#include "graphics/epsilon.h"

static void transform_face(std::vector<graphics::vec3>& points, graphics::ivec3& face, graphics::vec3* Mwt);

static void read_geom_image(
    const std::string& fname, std::vector<graphics::vec3>& points, 
    std::vector<graphics::vec3>& norms, 
    int& w, int& h, int& step)
{
	points.clear();
	norms.clear();
	std::string pos_file_name = fname + ".dat";
	std::string normal_file_name = fname + ".nor";

	FILE* fp = fopen(pos_file_name.c_str(), "rb");
	fread(&w, sizeof(int), 1, fp);
	fread(&h, sizeof(int), 1, fp);
	fread(&step, sizeof(int), 1, fp);
	int cnt = w * h;
	for (int i = 0; i < cnt; i++) {
		float v1, v2, v3;
		fread(&v1, sizeof(float), 1, fp);
		fread(&v2, sizeof(float), 1, fp);
		fread(&v3, sizeof(float), 1, fp);
		graphics::vec3 p(v1, v2, v3);
		points.push_back(p);
	}
	fclose(fp);
	fp = fopen(normal_file_name.c_str(), "rb");
	fread(&w, sizeof(int), 1, fp);
	fread(&h, sizeof(int), 1, fp);
	fread(&step, sizeof(int), 1, fp);
	cnt = w * h;
	for (int i = 0; i < cnt; i++) {
		float v1, v2, v3;
		fread(&v1, sizeof(float), 1, fp);
		fread(&v2, sizeof(float), 1, fp);
		fread(&v3, sizeof(float), 1, fp);
		graphics::vec3 p(v1, v2, v3);
		norms.push_back(p);
	}
	fclose(fp);
}

static void generate_face(
    const std::vector<graphics::vec3>& points,
    std::vector<graphics::ivec3>& faces,
    int& w, int& h)
{
    //(u, v) : generate 2 faces.
    // -> 0 :: (u, v-1) / 1 :: (u, v) / 2 :: (u+1, v-1) 
    // -> - :: (u+1, v) / 1 :: (u+1, v-1) / 2 :: (u, v)  
    int cnt = w * h;
    for(int v = 0; v < h -1 ; v++){
        for(int u = 1; u < w; u++){
            int p = v * w + u;

            if ( (points[p-1] > -1000.0 && points[p-1] > -1000.0 && points[p-1] > -1000.0) &&
            (points[p] > -1000.0 && points[p] > -1000.0 && points[p] > -1000.0) &&
            (points[p+w-1] > -1000.0 && points[p+w-1] > -1000.0 && points[p+w-1] > -1000.0))
            {
                faces.push_back(graphics::ivec3(p - 1, p , p + w - 1) );
            }

            if ( (points[p+w] > -1000.0 && points[p+w] > -1000.0 && points[p+w] > -1000.0) &&
            (points[p+w-1] > -1000.0 && points[p+w-1] > -1000.0 && points[p+w-1] > -1000.0) &&
            (points[p] > -1000.0 && points[p] > -1000.0 && points[p] > -1000.0))
            {        
                faces.push_back(graphics::ivec3(p + w, p + w - 1, p) );
            }
            //if((p + w) >= 2048 * 2048 || (p - 1) < 0) printf("error %d %d \n", v, u);
        }
    }

    //printf("w: %d h : %d\n", w, h);
    //printf("%d %d %d %d\n", faces.size(), faces[0][0], faces[0][1], faces[0][2]);

}

static void Cross( graphics::vec3& v1, graphics::vec3& v2, graphics::vec3& result)
{
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v2[0] * v1[2] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

static void normalize( graphics::vec3& v)
{
    float norm = sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2])) + 1e-8;
    v[0] /= norm;
    v[1] /= norm;
    v[2] /= norm;
}

static void transform_tangent_normal(
    std::vector<graphics::vec3>& points, 
    std::vector<graphics::vec3>& norms,
    std::vector<graphics::ivec3>& faces, 
    int& w, int& h)
{

    graphics::vec3 *tangent_norm = (graphics::vec3*) calloc(sizeof(graphics::vec3), norms.size());
    
    // initialize 3x3 matrix.
    graphics::vec3 Mwt[3];

    
    int len = faces.size();
    for (int i = 0; i < 3/*len*/ ; i++)
    {
        graphics::ivec3 f = faces[i];
        transform_face(points, f, Mwt);
    
        for(int ver = 0; ver < 3; ver++)
        {
            int ver_idx = f[ver];
            graphics::vec3 nxyz = norms[ver_idx];
            // [nx, ny, nz] * Mwt(3 x 3) = [tangent_nx, tangent_ny, tangent_nz]
            tangent_norm[ver_idx][0] += nxyz[0] * Mwt[0][0] + nxyz[1] * Mwt[1][0] + nxyz[2] * Mwt[2][0];
            tangent_norm[ver_idx][1] += nxyz[0] * Mwt[0][1] + nxyz[1] * Mwt[1][1] + nxyz[2] * Mwt[2][1];
            tangent_norm[ver_idx][2] += nxyz[0] * Mwt[0][2] + nxyz[1] * Mwt[1][2] + nxyz[2] * Mwt[2][2];
        }
        
    }

    for(int i=0; i< norms.size(); i++)
    {
        normalize(tangent_norm[i]);
        norms[i] = tangent_norm[i];
    }
    free(tangent_norm);
    
}



static void transform_face(std::vector<graphics::vec3> &points, graphics::ivec3 &face, graphics::vec3* Mwt)
{
    graphics::vec3 p0 = points[face[0]];
    graphics::vec3 p1 = points[face[1]];
    graphics::vec3 p2 = points[face[2]];

    float v1x = p1[0] - p0[0];
    float v1y = p1[1] - p0[1];
    float v1z = p1[2] - p0[2];
 
    float v2x = p2[0] - p0[0];
    float v2y = p2[1] - p0[1];
    float v2z = p2[2] - p0[2];

    float u1x = (v1x >= 0) ? 1 : -1;
    float u1y = 0;

    float u2x = 0;
    float u2y = (v1x >= 0) ? -1 : 1;

    // M :: t (M[0]), b (M[1]), n (M[2]) of surface (3 row vector)

    // Tangent :: (v1 * uv2y - v2 * uv1y) / det
    // Bitangent :: (v2 * uv1x - v1 * uv2x) / det

    // Tangent and Bitangent is not perpendicular !!
    Mwt[0][0] = (v1x * u2y - v2x * u1y);
    Mwt[0][1] = (v1y * u2y - v2y * u1y);
    Mwt[0][2] = (v1z * u2y - v2z * u1y);

    Mwt[1][0] = (v2x * u1x - v1x * u2x);
    Mwt[1][1] = (v2y * u1x - v1y * u2x);
    Mwt[1][2] = (v2z * u1x - v1z * u2x);

    // n : t cross b
    Cross(Mwt[0], Mwt[1], Mwt[2]);
    
    // normalize each vector
    normalize(Mwt[0]);
    normalize(Mwt[1]);
    normalize(Mwt[2]);

    // make b perpendicular.
    float dot = Mwt[0][0] * Mwt[1][0] + Mwt[0][1] * Mwt[1][1] + Mwt[0][2] * Mwt[1][2]; // t, b
    // b - t * dot (t, b)
    Mwt[1][0] -= dot * Mwt[0][0];
    Mwt[1][1] -= dot * Mwt[0][1];
    Mwt[1][2] -= dot * Mwt[0][2];
    normalize(Mwt[1]);

    // Headness of n (cross(n,t) dot b >= 0)
    graphics::vec3 n_t;
    Cross(Mwt[2], Mwt[0], n_t);
    float dot_n_t = n_t[0] * Mwt[1][0] + n_t[1] * Mwt[1][1] + n_t[2] * Mwt[1][2];
    if(dot_n_t < 0.0f)
    {
            Mwt[0][0] = -Mwt[0][0]; Mwt[0][1] = -Mwt[0][1]; Mwt[0][2] = -Mwt[0][2];
    }

/* test 
    printf("%f %f %f\n",p0[0], p0[1], p0[2]);
    printf("%f %f %f\n",p1[0], p1[1], p1[2]);
    printf("%f %f %f\n",p2[0], p2[1], p2[2]);
    printf("t : %f %f %f\nb : %f %f %f\nn : %f %f %f\n", 
    Mwt[0][0], Mwt[0][1], Mwt[0][2],
    Mwt[1][0], Mwt[1][1], Mwt[1][2],
    Mwt[2][0], Mwt[2][1], Mwt[2][2]);

    printf("dot test\n t, b\n%f\n t, n\n%f\nb, n\n%f\n",
    Mwt[0][0] * Mwt[1][0] + Mwt[0][1] * Mwt[1][1] + Mwt[0][2] * Mwt[1][2],
    Mwt[0][0] * Mwt[2][0] + Mwt[0][1] * Mwt[2][1] + Mwt[0][2] * Mwt[2][2],
    Mwt[1][0] * Mwt[2][0] + Mwt[1][1] * Mwt[2][1] + Mwt[1][2] * Mwt[2][2]);

*/
    // Transpose
    float M01 = Mwt[1][0], M02 = Mwt[2][0], M12 = Mwt[2][1];
    Mwt[1][0] = Mwt[0][1];
    Mwt[2][0] = Mwt[0][2];
    Mwt[2][1] = Mwt[1][2];
    Mwt[0][1] = M01;
    Mwt[0][2] = M02;
    Mwt[1][2] = M12;

}

static bool write_OBJfile(
    const char* out_obj, const char* in_texture, 
    std::vector<graphics::vec3>& points, 
    std::vector<graphics::vec3>& norms, 
    std::vector<graphics::ivec3>& faces,
    int w, int h, int precision
)
{
        FILE *f = fopen(out_obj, "wt");
        char formv[256];
        char formt[256];
        char formn[256];
        sprintf(formv, "v %%.%df %%.%df %%.%df\n",  precision, precision, precision);
        sprintf(formt, "vt %%.%df %%.%df\n",        precision, precision);
        sprintf(formn, "vn %%.%df %%.%df %%.%df\n", precision, precision, precision);
        
        if (f == NULL) 
        { 
            printf("Could not open file"); 
            return false; 
        }

        // fix mtl file name : texture.mtl 
        fprintf(f, "mtllib texture.mtl\n");

        std::vector<bool> mask;
        std::vector<int> mask_cnt; // maskcnt[i] points were masked out before vertex i.
        mask_cnt.push_back(0);
        for(int v=0; v<=points.size(); v++)
        {   
                if (points[v][0] > -1000.0 && points[v][1] > -1000.0 && points[v][2] > -1000.0) {
                    fprintf(f, formv, points[v][0], points[v][1], points[v][2]);
                    mask.push_back(true);
                    mask_cnt.push_back(mask_cnt.back());		
		        }
                else{
                    mask.push_back(false);
                    mask_cnt.push_back(mask_cnt.back() + 1);
                }
                                
        }
        int cnt = 0;
        for(int v=0; v < h; v++)
        {
            for(int u=0; u<w; u++)
            {
                if(mask[cnt]){
                    fprintf(f, formt, (float)u/w, 1-(float)v/h);
                }
                cnt++;
            }
        }
        cnt = 0;        
        for(int v=0; v<=norms.size(); v++)
        {       
                if(mask[cnt]){
                    fprintf(f, formn, norms[v][0], norms[v][1], norms[v][2]);
                }
                cnt++;
        }

        fprintf(f, "usemtl texture\n");
        for(int nF=0; nF<faces.size(); nF++)
        {       
                int i1 = faces[nF][0] + 1 - mask_cnt[faces[nF][0]]; 
                int i2 = faces[nF][1] + 1 - mask_cnt[faces[nF][1]]; 
                int i3 = faces[nF][2] + 1 - mask_cnt[faces[nF][2]];
                fprintf(f, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
                i1, i1, i1,
                i2, i2, i2,
                i3, i3, i3);   
        } 

        fclose(f);

        // Write Mtl file. (texture.mtl)
        FILE *ff = fopen("texture.mtl", "wt");
        if(ff == NULL) return false;
        fprintf(ff, "newmtl texture\n");
        fprintf(ff, "map_Kd %s\n", in_texture);
        fclose(ff);

        return true;
}


static bool write_PLYfile(
    const char* out_ply, 
    std::vector<graphics::vec3>& points, 
    std::vector<graphics::vec3>& points_normals, 
    std::vector<graphics::vec3>& image_points_color,
    std::vector<graphics::ivec3>& faces,
    int filetype)
{
	FILE *fp;
	fp = fopen(out_ply, "wb");
	if (fp == NULL) {
		//LOG("fail to open ply file \n");
		return false;
	}
	fprintf(fp, "ply\n");
	if (filetype == 0)
		fprintf(fp, "format ascii 1.0\n");
	else
    {
        fprintf(fp, "format binary_little_endian 1.0\n");
    }
    
	fprintf(fp, "element vertex %d\n", points.size());
	fprintf(fp, "property float x \n");
	fprintf(fp, "property float y \n");
	fprintf(fp, "property float z \n");
	fprintf(fp, "property float nx \n");
	fprintf(fp, "property float ny \n");
	fprintf(fp, "property float nz \n");
	fprintf(fp, "property uchar red \n");
	fprintf(fp, "property uchar green \n");
	fprintf(fp, "property uchar blue \n");
	fprintf(fp, "element face 0 \n");
	fprintf(fp, "property list uchar int vertex_indices \n");
	fprintf(fp, "end_header\n");
	float pt[3], n[3];
	unsigned char c[3];
	void* data;
	for (int i = 0; i < points.size(); i++)
	{
		pt[0] = points[i][0];
		pt[1] = points[i][1];
		pt[2] = points[i][2];
		if (_isnan(points_normals[i][0]) || _isnan(points_normals[i][1]) || _isnan(points_normals[i][2]))
		{
			points_normals[i] = graphics::vec3(0, 0, 1);
		}
        /*
		if (apx_equal(norm(points_normals[i]), 0, zero_epsilon)) {
			points_normals[i] = vec3(0, 0, 1);
		}*/
		n[0] = points_normals[i][0];
		n[1] = points_normals[i][1];
		n[2] = points_normals[i][2];
		//n[0] = 0.0; n[1] = 0.0; n[2] = 0.0;
		c[0] = image_points_color[i][0] * 255;
		c[1] = image_points_color[i][1] * 255;
		c[2] = image_points_color[i][2] * 255;
		if (filetype == 0)
			fprintf(fp, "%E %E %E %E %E %E %u %u %u\n", pt[0], pt[1], pt[2], n[0], n[1], n[2], c[0], c[1], c[2]);
		else {
			data = pt;
			fwrite(data, sizeof(float), 3, fp);
			data = n;
			fwrite(data, sizeof(float), 3, fp);
			data = c;
			fwrite(data, sizeof(unsigned char), 3, fp);
		}
	}
    unsigned char three = 3;
    int f[3];
    for( int i = 0; i < faces.size(); i++)
    {
        f[0] = faces[i][0];
        f[1] = faces[i][1];
        f[2] = faces[i][2];
        if(filetype == 0)
        {
            fprintf(fp, "%u %d %d %d\n", 3, f[0], f[1], f[2]);
        }
        else
        {
            data = &three;
            fwrite(data, sizeof(unsigned char), 1, fp);
            data = f;
            fwrite(data, sizeof(int), 3, fp);
        }
        

    }
	fclose(fp);
	return true;
}

int main()
{
    std::string infile = "./bin/geom_img";
    std::vector<graphics::vec3> points;
    std::vector<graphics::vec3> norms;
    std::vector<graphics::ivec3> faces;
    std::vector<graphics::vec3> color;

    int w, h, step;
    int precision = 6;

    read_geom_image(infile, points, norms, w, h, step);
    generate_face(points, faces, w, h);
    //transform_tangent_normal(points, norms, faces, w, h);
    std::string output = "out.obj";
    std::string input_texture = "blended_texture.bmp";

    write_OBJfile(output.c_str(), input_texture.c_str(), points, norms, faces, w, h, precision);
    //write_PLYfile(output.c_str(), points, norms, color, faces, 1);
}   