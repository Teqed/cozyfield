use bevy::{
    math::vec2,
    prelude::*,
    render::{
        settings::{Backends, RenderCreation, WgpuSettings},
        RenderPlugin,
    },
    sprite::{MaterialMesh2dBundle, Mesh2dHandle},
};
#[derive(Component)]
struct Velocity {
    speed_x: f32,
    speed_y: f32,
}
#[derive(Component)]
struct GridVelocity;
#[derive(Clone)]
struct SubGrid(Vec<Vec2>);
#[derive(Clone)]
struct Grid(Vec<SubGrid>);
impl std::ops::Index<usize> for SubGrid {
    type Output = Vec2;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}
impl std::ops::IndexMut<usize> for SubGrid {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}
impl std::ops::Index<usize> for Grid {
    type Output = SubGrid;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}
impl std::ops::IndexMut<usize> for Grid {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}
const X_EXTENT: f32 = 500.;
const Y_EXTENT: f32 = 500.;
fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
) {
    let mut camera = Camera2dBundle::default();
    camera.transform.translation.x = 0.;
    commands.spawn(camera);
    let circle = Mesh2dHandle(meshes.add(Circle { radius: 25.0 }));
    let circle2 = Mesh2dHandle(meshes.add(Circle { radius: 25.0 }));
    let circle3 = Mesh2dHandle(meshes.add(Circle { radius: 25.0 }));
    let transform = Transform::from_xyz(0., 0., 0.);
    let transform2 = Transform::from_xyz(100., 100., 1.);
    let transform3 = Transform::from_xyz(0., 200., 2.);
    let color = Color::rgb(1.0, 0.0, 0.0);
    let color2 = Color::rgb(0.0, 0.0, 1.0);
    let color3 = Color::rgb(0.0, 1.0, 0.0);
    commands.spawn((
        MaterialMesh2dBundle {
            mesh: circle,
            material: materials.add(color),
            transform,
            ..Default::default()
        },
        Velocity {
            speed_x: 0.1,
            speed_y: 0.,
        },
    ));
    commands.spawn((
        MaterialMesh2dBundle {
            mesh: circle2,
            material: materials.add(color2),
            transform: transform2,
            ..Default::default()
        },
        Velocity {
            speed_x: 0.01,
            speed_y: -0.2,
        },
    ));
    commands.spawn((
        MaterialMesh2dBundle {
            mesh: circle3,
            material: materials.add(color3),
            transform: transform3,
            ..Default::default()
        },
        Velocity {
            speed_x: 0.,
            speed_y: 0.1,
        },
    ));
    // Spawn 20 randomcly placed and colored dots
    for _ in 0..20 {
        let mut dot_transform = Transform::from_xyz(
            (rand::random::<f32>() - 0.5) * X_EXTENT,
            (rand::random::<f32>() - 0.5) * Y_EXTENT,
            0.,
        );
        dot_transform.scale = Vec3::splat(rand::random::<f32>());
        let dot = Mesh2dHandle(meshes.add(Circle { radius: 25. }));
        let hue = rand::random::<f32>() * 360.;
        let dot_color = Color::hsl(hue, 1., 0.5);
        commands.spawn((
            MaterialMesh2dBundle {
                mesh: dot,
                material: materials.add(dot_color),
                transform: dot_transform,
                ..Default::default()
            },
            Velocity {
                speed_x: 0.,
                speed_y: 0.,
            },
        ));
    }
    commands.spawn(System::new());
}
#[derive(Component)]
struct System {
    peaks: Vec<Peak>,               // Positions of the centers of divergence
                               // and their amplitudes
    grid: Grid,                // 3D array to store the combined Gaussian influences
    resolution: Resolution,    // Defines the resolution of the grid (width, height)
    cell_size: f32,            // The physical size of each cell in the grid
}
struct Peak {
    grid_x: usize,
    grid_y: usize,
    amplitude: f32,
}
struct Resolution {
    width: usize,
    height: usize,
}
const DIMENSION : usize = 100;
const CELL : f32 = 10.;
impl System {
    fn new() -> Self {
        Self {
            peaks: Vec::new(),
            grid: Grid(vec![SubGrid(vec![Vec2 { x: 0., y: 0. }; DIMENSION]); DIMENSION]),
            resolution: Resolution {
                width: DIMENSION,
                height: DIMENSION,
            },
            cell_size: CELL,
        }
    }
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_precision_loss,
        clippy::cast_sign_loss
    )]
    #[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss, clippy::cast_sign_loss)]
    fn position_to_cell(&self, position: f32) -> usize {
        if position < -X_EXTENT {
            0
        } else if position > X_EXTENT {
            self.resolution.width - 1
        } else {
            ((position + X_EXTENT) / self.cell_size) as usize
        }
    }
    fn collect_peak(&mut self, peak_x: f32, peak_y: f32, peak_amplitude: f32) {
        self.peaks.push(Peak {
            grid_x: self.position_to_cell(peak_x),
            grid_y: self.position_to_cell(peak_y),
            amplitude: peak_amplitude,
        });
    }
    fn clear_peaks(&mut self) {
        self.peaks.clear();
        // self.grid = Grid(vec![SubGrid(vec![Vec2 { x: 0., y: 0. }; 100]); 100]);        
        // Have the boundaries absorb the heat
        // The heat equation is a partial differential equation which describes the flow of heat in a given domain. For a 2D domain, it is given by the following equation:
        // ∂u/∂t = α(∂²u/∂x² + ∂²u/∂y²)
        // where u is the temperature, t is the time, x and y are the spatial coordinates, and α is the thermal diffusivity.
        // The equation is solved numerically using the finite difference method, which consists in discretizing the domain into a grid and approximating the derivatives by finite differences. The discretization scheme used here is the Crank-Nicolson method, which is unconditionally stable and second-order accurate in time and space.
        // The resulting system of equations is then solved using the Conjugate Gradient method.
        for _ in 0..10 {
            for x in 0..self.resolution.width {
                self.grid[0][x] = Vec2::new(-0.1, -0.1);
                self.grid[self.resolution.height - 1][x] = Vec2::new(-0.1, -0.1);
            }
            for y in 0..self.resolution.height {
                self.grid[y][0] = Vec2::new(-0.1, -0.1);
                self.grid[y][self.resolution.width - 1] = Vec2::new(-0.1, -0.1);
            }
            let thermal_diffusivity = 0.1; // Thermal diffusivity
            let dt = 2.;   // Time step
            let dx = 1_f32; // Spatial step
            let dy = 1_f32; // Spatial step
            let alpha = thermal_diffusivity * dt / (dx * dy);
            let mut new_grid = self.grid.clone();
            for x in 1..self.resolution.width - 1 {
                for y in 1..self.resolution.height - 1 {
                    let laplacian = (self.grid[y][x + 1] + self.grid[y][x - 1] - 2. * self.grid[y][x])
                        / dx.powi(2)
                        + (self.grid[y + 1][x] + self.grid[y - 1][x] - 2. * self.grid[y][x])
                            / dy.powi(2);
                    new_grid[y][x] = self.grid[y][x] + alpha * laplacian;
                }
            }
            self.grid = new_grid;
        }
        // Cool down the grid
        for x in 0..self.resolution.width {
            for y in 0..self.resolution.height {
                if self.grid[y][x].length() > 0. {
                    self.grid[y][x] *= 0.99;
                }
            }
        }
    }
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_precision_loss,
        clippy::cast_sign_loss
    )]
    fn compute_influence(&mut self) {
        for peak in &self.peaks {
            // Find the grid cell that contains the peak
            let amplitude = peak.amplitude;
            // Deposit the peak's influence into the grid cell
            self.grid[peak.grid_y][peak.grid_x] += vec2(amplitude, amplitude);
        }
    }
}
#[allow(clippy::needless_pass_by_value)]
fn move_shapes(
    mut query: Query<(&mut Transform, &mut Velocity, &Mesh2dHandle), With<Velocity>>,
    system: Query<&mut System>,
) {
    for mut entity in &mut query {
        let mut delta = vec2(0., 0.);
        let cell_x = system.single().position_to_cell(entity.0.translation.x);
        let cell_y = system.single().position_to_cell(entity.0.translation.y);
        if cell_x < 1 || cell_x >= system.single().resolution.width - 1
            || cell_y < 1
            || cell_y >= system.single().resolution.height - 1
        {
            // If the shape is outside the grid, move it back to the center
            entity.1.speed_x -= entity.0.translation.x * 0.001;
            entity.1.speed_y -= entity.0.translation.y * 0.001;
            entity.0.translation.x += entity.1.speed_x;
            entity.0.translation.y += entity.1.speed_y;
            continue;
        }
        let left_influence = system.single().grid[cell_y][cell_x - 1];
        let right_influence = system.single().grid[cell_y][cell_x + 1];
        let up_influence = system.single().grid[cell_y - 1][cell_x];
        let down_influence = system.single().grid[cell_y + 1][cell_x];
        delta.x = (right_influence - left_influence).x;
        delta.y = (down_influence - up_influence).y;
        let mass = entity.0.scale.x;
        entity.1.speed_x += delta.x / mass;
        entity.1.speed_y += delta.y / mass;
        // Apply a slight damping to the speed
        entity.1.speed_x *= 0.99999;
        entity.1.speed_y *= 0.99999;
        entity.0.translation.x += entity.1.speed_x;
        entity.0.translation.y += entity.1.speed_y;
        if entity.0.translation.x > X_EXTENT - 10. || entity.0.translation.x < -X_EXTENT + 10. {
            entity.1.speed_x = -entity.1.speed_x * 0.8;
        }
        if entity.0.translation.y > Y_EXTENT - 10. || entity.0.translation.y < -Y_EXTENT + 10. {
            entity.1.speed_y = -entity.1.speed_y * 0.8;
        }
        if entity.0.translation.x >= X_EXTENT {
            entity.0.translation.x = X_EXTENT - 1.;
        } else if entity.0.translation.x <= -X_EXTENT {
            entity.0.translation.x = -X_EXTENT + 1.;
        }
        if entity.0.translation.y >= Y_EXTENT {
            entity.0.translation.y = Y_EXTENT - 1.;
        } else if entity.0.translation.y <= -Y_EXTENT {
            entity.0.translation.y = -Y_EXTENT + 1.;
        }
    }
}
#[allow(clippy::needless_pass_by_value)]
fn add_peaks(
    mut system: Query<&mut System>,
    query: Query<(&Transform, &Mesh2dHandle), With<Velocity>>,
) {
    system.single_mut().clear_peaks();
    for entity in &mut query.iter() {
        let radius = entity.0.scale.x;
        system
            .single_mut()
            .collect_peak(entity.0.translation.x, entity.0.translation.y, radius);
    }
    system.single_mut().compute_influence();
}
#[allow(clippy::needless_pass_by_value, clippy::cast_precision_loss)]
fn change_dot_color(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    system: Query<&System>,
    dots: Query<(Entity, &GridVelocity)>,
) {
    // Remove all old dots
    for (entity, _) in &mut dots.iter() {
        commands.entity(entity).despawn();
    }
    // for x in -10..10 {
    //     for y in -10..10 {
    //         let translation = Vec3::new(x as f32 * 50., y as f32 * 50., -10.);
    for x in 0..system.single().resolution.width {
        if x % 4 != 0 {
            continue;
        }
        for y in 0..system.single().resolution.height {
            if y % 4 != 0 {
                continue;
            }
            let translation = Vec3::new(
                (x as f32 - system.single().resolution.width as f32 / 2.).mul_add(system.single().cell_size, system.single().cell_size / 2.),
                (y as f32 - system.single().resolution.height as f32 / 2.).mul_add(system.single().cell_size, system.single().cell_size / 2.),
                -10.,
            );
            let grid_velocity = system.single().grid[y][x];
            let radius = grid_velocity.length().mul_add(10., 5.);
            let dot_radius = radius.abs().clamp(0.2, 20.);
            let dot = Mesh2dHandle(meshes.add(Circle { radius: dot_radius }));
            let hue = (grid_velocity.y.atan2(grid_velocity.x) + std::f32::consts::PI)
                / (2. * std::f32::consts::PI);
            let hue = hue * 360.;
            let saturation = // out of 1.0
                1. - (grid_velocity.length() / 100.).min(0.01);
            let color = Color::hsl(hue, saturation, saturation * 0.3);
            commands.spawn((
                MaterialMesh2dBundle {
                    mesh: dot,
                    material: materials.add(color),
                    transform: Transform::from_translation(translation),
                    ..Default::default()
                },
                GridVelocity {},
            ));
        }
    }
}
struct PluginBackendProvided;
impl Plugin for PluginBackendProvided {
    fn build(&self, app: &mut App) {
        app.add_plugins(DefaultPlugins.set(RenderPlugin {
            render_creation: RenderCreation::Automatic(WgpuSettings {
                backends: Some(Backends::VULKAN),
                ..default()
            }),
            synchronous_pipeline_compilation: false,
        }));
    }
}
fn main() {
    App::new()
        .add_plugins(PluginBackendProvided)
        .add_systems(Startup, setup)
        .add_systems(Update, (add_peaks, move_shapes.after(add_peaks)))
        .add_systems(Update, change_dot_color.after(add_peaks))
        .run();
}
