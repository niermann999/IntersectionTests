
namespace vec_intr {

namespace g_bench {

//----------------------------------------------------Fill Data Types

// Vertical data as in detray
void VertSetup::SetUp(const ::benchmark::State& state) {
      // Only the first thread sets the test env up
      if (state.thread_index != 0) return;

      ray_dir.obj   = vector_s::type::Random();
      ray_point.obj = vector_s::type::Random();

      planes.reserve(state.range(0));
      for (size_t i = 0; i < state.range(0); i++) {
        pl_normals = {.obj = vector_s::type::Random()};
        pl_points  = {.obj = vector_s::type::Random()};

        planes.push_back({.normals = pl_normals, .points = pl_points});
      }

      // vertical data containers
      ray = {.direction = std::move(ray_dir), 
             .point     = std::move(ray_point)};
    }

    void VertSetup::TearDown(const ::benchmark::State& state) {
      // Only one thread frees recources
      if (state.thread_index == 0) {
        planes.clear();
      }
    }


// Vertical data contained in structs
    void HybridSetup::SetUp(const ::benchmark::State& state) {
      // Only the first thread sets the test env up
      if (state.thread_index != 0) return;

      // AoS data
      pl_normals_struct.reserve(state.range(0));
      pl_points_struct.reserve(state.range(0));
      for (size_t i = 0; i < state.range(0); i++) {
        pl_normals_struct.push_back({.x=uni(), .y=uni(), .z=uni()});
        pl_points_struct.push_back( {.x=uni(), .y=uni(), .z=uni()});
      }

      // Horizontal data (interleaved)
      ray_dir_hor.obj   = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
      ray_point_hor.obj = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};

      // AoS container
      ray_struct    = {.direction = std::move(ray_dir_hor), 
                       .point     = std::move(ray_point_hor)}; 
      planes_struct = {.normals = std::move(pl_normals_struct), 
                       .points  = std::move(pl_points_struct)};
    }

    void HybridSetup::TearDown(const ::benchmark::State& state) {
      // Only one thread frees recources
      if (state.thread_index == 0) {
        planes_struct.normals.clear();
        planes_struct.points.clear();
      }
    }


// Horizontal data as interleaved horizonral vectors
    void HorizSetup::SetUp(const ::benchmark::State& state) {
      // Only the first thread sets the test env up
      if (state.thread_index != 0) return;

      // Horizontal data (interleaved)
      ray_dir_hor.obj   = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
      ray_point_hor.obj = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};

      // horizontal ray container
      ray_hor = {.direction = std::move(ray_dir_hor),    
                 .point     = std::move(ray_point_hor)};

      planes_hor.reserve(6* state.range(0)/Scalar_v::Size + 1);
      for (size_t s = 0; s < state.range(0)/Scalar_v::Size; s++) {
        pl_normals_hor.obj = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
        pl_points_hor.obj  = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
        planes_hor.push_back({.normals = std::move(pl_normals_hor), 
                              .points  = std::move(pl_points_hor)});
      }
      // padding at the end of data container needed (just add one more calculation for simplicity)
      if (state.range(0)/Scalar_v::Size != 0) {
        pl_normals_hor.obj = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
        pl_points_hor.obj  = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
        planes_hor.push_back({.normals = std::move(pl_normals_hor), 
                              .points  = std::move(pl_points_hor)});
      }
    }

    void HorizSetup::TearDown(const ::benchmark::State& state) {
      // Only one thread frees recources
      if (state.thread_index == 0) {
        planes_hor.clear();
      }
    }

} // namespace g_bench

} //namespace vec_intr