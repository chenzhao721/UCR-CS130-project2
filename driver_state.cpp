#include "driver_state.h"
#include <math.h>
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    //state.image_color=0;
    //state.image_depth=0;
    
    //my code here
    state.image_color  = (pixel *)malloc(width*height*sizeof(pixel));
    state.image_depth  = (float *)malloc(width*height*sizeof(float));
    //state.image_color = new pixel[state.image_width][state.image_height];
    
    for (int i = 0; i < width*height; i++){
        state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = 20000.0;
    }
    
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    //my code here
    data_geometry *array = (data_geometry *)malloc(state.num_vertices*sizeof(data_geometry));
    for (int i = 0; i < state.num_vertices; i++)  array[i].data = (float *)malloc(state.floats_per_vertex*sizeof(float));
    
    data_vertex in;
    for (int i = 0; i < state.num_vertices; i++){ 
        in.data = &state.vertex_data[i*state.floats_per_vertex];
        (*state.vertex_shader)(in, array[i], state.uniform_data);
        array[i].data = in.data; 
        
    }
    
    switch(type)
    {
        case render_type::triangle: 
            //std::cout<<state.num_triangles<<std::endl;
            for (int i = 0; i < floor(state.num_vertices/3.0); i++){
                clip_triangle(state, array[3*i], array[3*i+1], array[3*i+2],0); 
                
            }
            break;
        case render_type::indexed:  
            for (int i = 0; i < state.num_triangles; i++){
                clip_triangle(state, array[state.index_data[3*i]], array[state.index_data[3*i+1]], array[state.index_data[3*i+2]],0); 
            }
            break;
        case render_type::fan: 
            for (int i = 1; i <= (state.num_vertices-2); i++){
                clip_triangle(state, array[0], array[i], array[i+1],0); 
            }
            break;
            
        case render_type::strip:
            for (int i = 2; i <= state.num_vertices-1; i++){
                if(i%2==1)
                    clip_triangle(state, array[i-2], array[i-1], array[i],0); 
                else
                    clip_triangle(state, array[i-1], array[i-2], array[i],0); 
            }
            break;
    }
           
 
        
    //std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        //std::cout<<"clipping finished,face:"<<face<<std::endl;
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    
    float alpha,beta;
    float a,b,c;
    data_geometry ab,ac,bc;
    float wa,wb,wc;
    bool flag = 1;
    
    //std::cout<<"v0.0 :"<<v0.gl_Position[0]<<" v0.1 :"<<v0.gl_Position[1]<<" v0.2 :"<<v0.gl_Position[2]<<std::endl;
    a = v0.gl_Position[face/2] ;
    b = v1.gl_Position[face/2] ;
    c = v2.gl_Position[face/2] ;
    wa = v0.gl_Position[3] ;
    wb = v1.gl_Position[3] ;
    wc = v2.gl_Position[3] ;
    
    ab.data = new float[state.floats_per_vertex];
    ac.data = new float[state.floats_per_vertex];
    bc.data = new float[state.floats_per_vertex];
    
    
    if(face==1 or face==3 or face==5){
        wa = -wa;
        wb = -wb;
        wc = -wc;
        flag = 0;
    }
    
    if((a<=wa and b<=wb and c<=wc and (flag == 1)) or (a>=wa and b>=wb and c>=wc and (flag == 0))){
        //std::cout<<"no clippng,face:"<<face<<std::endl;
        clip_triangle(state, v0, v1, v2,face+1);
        
        return;
        
    }
    
    if((a>wa and b>wb and c>wc and (flag == 1)) or (a<wa and b<wb and c<wc and (flag == 0))){
        
        return;
    }
    
    if((a<=wa and b>wb and c>wc and (flag == 1)) or (a>=wa and b<wb and c<wc and (flag == 0))){        //a in, bc out
       
       
       alpha = (wb-b)/(a-wa+wb-b);
       beta= (wc-c)/(a-wa+wc-c);
       
       for(int i = 0; i<state.floats_per_vertex; i++){
           ab.data[i] = alpha*v0.data[i] + (1-alpha)*v1.data[i];
           ac.data[i] = beta*v0.data[i] + (1-beta)*v2.data[i];
           
           if(state.interp_rules[i]==interp_type::noperspective){
               float wp = alpha*wa+(1-alpha)*wb;
               float wp2 = beta*wa+(1-beta)*wc;
               float alpha_c = alpha*wa/wp;
               float beta_c = beta*wa/wp2;
               ab.data[i] = alpha_c*v0.data[i] + (1-alpha_c)*v1.data[i];
               ac.data[i] = beta_c*v0.data[i] + (1-beta_c)*v2.data[i];
           }
           
       } 
       for(int i = 0; i<4; i++){ 
           ab.gl_Position[i] = alpha*v0.gl_Position[i]  + (1-alpha)*v1.gl_Position[i] ;
           ac.gl_Position[i] = beta*v0.gl_Position[i]  + (1-beta)*v2.gl_Position[i] ;
       }
       
       clip_triangle(state, v0, ab, ac, face+1);
       return;
    }
    
    if((a<=wa and b<=wb and c>wc and (flag == 1)) or (a>=wa and b>=wb and c<wc and (flag == 0))){  //a,b in,  c out
       
       alpha = (wc-c)/(a-wa+wc-c);
       beta= (wc-c)/(b-wb+wc-c);     
       for(int i = 0; i<state.floats_per_vertex; i++){
           ac.data[i] = alpha*v0.data[i] + (1-alpha)*v2.data[i];
           bc.data[i] = beta*v1.data[i] + (1-beta)*v2.data[i];
           if(state.interp_rules[i]==interp_type::noperspective){
               float wp = alpha*wa+(1-alpha)*wc;
               float wp2 = beta*wb+(1-beta)*wc;
               float alpha_c = alpha*wa/wp;
               float beta_c = beta*wb/wp2;
               ac.data[i] = alpha_c*v0.data[i] + (1-alpha_c)*v2.data[i];
               bc.data[i] = beta_c*v1.data[i] + (1-beta_c)*v2.data[i];
           }
                      
           
       } 
       for(int i = 0; i<4; i++){ 
           ac.gl_Position[i] = alpha*v0.gl_Position[i]  + (1-alpha)*v2.gl_Position[i] ;
           bc.gl_Position[i] = beta*v1.gl_Position[i]  + (1-beta)*v2.gl_Position[i] ;
       }
       
       clip_triangle(state, v0, v1, ac, face+1);
       clip_triangle(state, ac, v1, bc, face+1);
       return;
    }
    
    if((a>wa and b<=wb and c>wc and (flag == 1)) or (a<wa and b>=wb and c<wc and (flag == 0))){    //b in, ac out
       
       alpha = (wa-a)/(b-wb+wa-a);
       beta= (wc-c)/(b-wb+wc-c);
       for(int i = 0; i<state.floats_per_vertex; i++){
           ab.data[i] = alpha*v1.data[i] + (1-alpha)*v0.data[i];
           bc.data[i] = beta*v1.data[i] + (1-beta)*v2.data[i];
           if(state.interp_rules[i]==interp_type::noperspective){
               float wp = alpha*wb+(1-alpha)*wa;
               float wp2 = beta*wb+(1-beta)*wc;
               float alpha_c = alpha*wb/wp;
               float beta_c = beta*wb/wp2;
               ab.data[i] = alpha_c*v1.data[i] + (1-alpha_c)*v0.data[i];
               bc.data[i] = beta_c*v1.data[i] + (1-beta_c)*v2.data[i];
           }
           
           
       } 
       for(int i = 0; i<4; i++){ 
           ab.gl_Position[i] = alpha*v1.gl_Position[i]  + (1-alpha)*v0.gl_Position[i] ;
           bc.gl_Position[i] = beta*v1.gl_Position[i]  + (1-beta)*v2.gl_Position[i] ;
       }
       
       clip_triangle(state, v1, bc, ab, face+1);
       return;
    }
    
    if((a>wa and b<=wb and c<=wc and (flag == 1)) or (a<wa and b>=wb and c>=wc and (flag == 0))){     //bc in, a out
       
       
       alpha = (wa-a)/(b-wb+wa-a);
       beta= (wa-a)/(c-wc+wa-a);
       
       for(int i = 0; i<state.floats_per_vertex; i++){
           
           ab.data[i] = alpha*v1.data[i] + (1-alpha)*v0.data[i];
           ac.data[i] = beta*v2.data[i] + (1-beta)*v0.data[i];
           if(state.interp_rules[i]==interp_type::noperspective){
               float wp = alpha*wb+(1-alpha)*wa;
               float wp2 = beta*wc+(1-beta)*wa;
               float alpha_c = alpha*wb/wp;
               float beta_c = beta*wc/wp2;
               ab.data[i] = alpha_c*v1.data[i] + (1-alpha_c)*v0.data[i];
               ac.data[i] = beta_c*v2.data[i] + (1-beta_c)*v0.data[i];
           }
           
           
       } 
       for(int i = 0; i<4; i++){ 
           
           ab.gl_Position[i] = alpha*v1.gl_Position[i]  + (1-alpha)*v0.gl_Position[i] ;
           ac.gl_Position[i] = beta*v2.gl_Position[i]  + (1-beta)*v0.gl_Position[i] ;
       }
       
       clip_triangle(state, v1, v2, ab, face+1);
       clip_triangle(state, ab, v2, ac, face+1);
       return;
    }
    
    if((a>wa and b>wb and c<=wc and (flag == 1)) or (a<wa and b<wb and c>=wc and (flag == 0))){   //c in, ab out
       
       
       alpha = (wa-a)/(c-wc+wa-a);
       beta= (wb-b)/(c-wc+wb-b);
       
       for(int i = 0; i<state.floats_per_vertex; i++){
           ac.data[i] = alpha*v2.data[i] + (1-alpha)*v0.data[i];
           bc.data[i] = beta*v2.data[i] + (1-beta)*v1.data[i];
           if(state.interp_rules[i]==interp_type::noperspective){
               float wp = alpha*wc+(1-alpha)*wa;
               float wp2 = beta*wc+(1-beta)*wb;
               float alpha_c = alpha*wc/wp;
               float beta_c = beta*wc/wp2;
               ac.data[i] = alpha_c*v2.data[i] + (1-alpha_c)*v0.data[i];
               bc.data[i] = beta_c*v2.data[i] + (1-beta_c)*v1.data[i];
           }
           
           
       } 
       for(int i = 0; i<4; i++){ 
           ac.gl_Position[i] = alpha*v2.gl_Position[i]  + (1-alpha)*v0.gl_Position[i] ;
           bc.gl_Position[i] = beta*v2.gl_Position[i]  + (1-beta)*v1.gl_Position[i] ;
       }
       
       clip_triangle(state, v2, ac, bc, face+1);
       return;
    }
    
    if((a<=wa and b>wb and c<=wc and (flag == 1)) or (a>=wa and b<wb and c>=wc and (flag == 0))){ //ac in, b out
       
       
       alpha = (wb-b)/(a-wa+wb-b);
       beta= (wb-b)/(c-wc+wb-b);
       
       for(int i = 0; i<state.floats_per_vertex; i++){
           ab.data[i] = alpha*v0.data[i] + (1-alpha)*v1.data[i];
           bc.data[i] = beta*v2.data[i] + (1-beta)*v1.data[i];
           if(state.interp_rules[i]==interp_type::noperspective){
               float wp = alpha*wa+(1-alpha)*wb;
               float wp2 = beta*wc+(1-beta)*wb;
               float alpha_c = alpha*wa/wp;
               float beta_c = beta*wc/wp2;
               ab.data[i] = alpha_c*v0.data[i] + (1-alpha_c)*v1.data[i];
               bc.data[i] = beta_c*v2.data[i] + (1-beta_c)*v1.data[i];
           }
           
           
       } 
       for(int i = 0; i<4; i++){ 
           ab.gl_Position[i] = alpha*v0.gl_Position[i]  + (1-alpha)*v1.gl_Position[i] ;
           bc.gl_Position[i] = beta*v2.gl_Position[i]  + (1-beta)*v1.gl_Position[i] ;
       }
       
       clip_triangle(state, v2, v0, bc, face+1);
       clip_triangle(state, bc, v0, ab, face+1);
       return;
    }
       
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    float Ai, Aj, Bi, Bj, Ci, Cj;
    float zf;                           //z-buffering temp
    int r,g,b;
    int width = state.image_width;
    int height = state.image_height;
    Ai = (width)*(v0.gl_Position[0]/v0.gl_Position[3] +1)/2.0-0.5;
    Aj = (height)*(v0.gl_Position[1]/v0.gl_Position[3] +1)/2.0-0.5;
    Bi = (width)*(v1.gl_Position[0]/v1.gl_Position[3] +1)/2.0-0.5;
    Bj = (height)*(v1.gl_Position[1]/v1.gl_Position[3] +1)/2.0-0.5;
    Ci = (width)*(v2.gl_Position[0]/v2.gl_Position[3] +1)/2.0-0.5;
    Cj = (height)*(v2.gl_Position[1]/v2.gl_Position[3] +1)/2.0-0.5;
    
    //int area_ABC = 0.5((Bi*Cj − Ci* Bj ) + (Ci* Aj − Ai* Cj ) + (Ai* Bj − Bi* Aj ));
    float area_ABC = (Bi-Ai)*(Cj-Aj)-(Bj-Aj)*(Ci-Ai);
    
    float alpha, beta, gamma;
    data_output color;
    data_fragment in;
    in.data = new float[state.floats_per_vertex];
            
    
    for (int i = 0; i < width; i++){
        for (int j = 0; j < height; j++){
            alpha = ((Bi-i)*(Cj-j)-(Bj-j)*(Ci-i))/area_ABC;
            beta = ((Ai-i)*(Bj-j)-(Aj-j)*(Bi-i))/area_ABC;
            gamma = ((Ci-i)*(Aj-j)-(Cj-j)*(Ai-i))/area_ABC;
            double wp = 1.0/(alpha/v0.gl_Position[3]+beta/v2.gl_Position[3]+gamma/v1.gl_Position[3]);
            float alpha_c = alpha*wp/v0.gl_Position[3];
            float beta_c = beta*wp/v2.gl_Position[3];
            float gamma_c = gamma*wp/v1.gl_Position[3];
            
            if(alpha>=0&&beta>=0&&gamma>=0){
                for(int i = 0; i<state.floats_per_vertex; i++){           
                    switch(state.interp_rules[i]){
                        case interp_type::flat:                      
                            in.data[i] = v0.data[i];
                        break;
                        
                        case interp_type::noperspective:                       
                            in.data[i] = alpha*v0.data[i]+beta*v2.data[i]+gamma*v1.data[i];     
                            //std::cout<<"noperspective, i is :"<<i<<std::endl;               
                        break;
                        
                        case interp_type::smooth:                             
                            //std::cout<<"alpha:"<<alpha<<" beta:"<<beta<<" gamma:"<<gamma<<std::endl;
                            //std::cout<<"sum:"<<(alpha+beta+gamma)<<std::endl;
                            in.data[i] = alpha_c*v0.data[i] + beta_c*v2.data[i] + gamma_c*v1.data[i]; 
                            //std::cout<<"smooth, i is :"<<i<<std::endl;           
                        break;
                    }
                }
                (*state.fragment_shader)(in, color, state.uniform_data);
                r = int(color.output_color[0]*255);
                g = int(color.output_color[1]*255);
                b = int(color.output_color[2]*255);
                
                zf = alpha*v0.gl_Position[2]/v0.gl_Position[3] + beta*v2.gl_Position[2]/v2.gl_Position[3] + gamma*v1.gl_Position[2]/v1.gl_Position[3];  //z-buffering
                //std::cout<<"zf is :"<<zf<<"|| ";
                
                if(zf < state.image_depth[j*width+i]){
                    state.image_color[j*width+i] = make_pixel(r, g, b);
                    state.image_depth[j*width+i] = zf;
                }
                
                //state.image_color[j*width+i] = make_pixel(r, g, b);
                
            }
            
            
        }
    }  
    
    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

