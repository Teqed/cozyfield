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
    peaks: Vec<Vec2>,          // Positions of the centers of divergence
    peak_amplitudes: Vec<f32>, // 'a' in the Gaussian function
    width: f32,                // 'c' in the Gaussian function
    grid: Grid,                // 3D array to store the combined Gaussian influences
    resolution: Resolution,    // Defines the resolution of the grid (width, height)
    cell_size: f32,            // The physical size of each cell in the grid
}
struct Resolution {
    width: usize,
    height: usize,
}
impl System {
    fn new() -> Self {
        Self {
            peaks: Vec::new(),
            peak_amplitudes: Vec::new(),
            width: 75.,
            grid: Grid(vec![SubGrid(vec![Vec2 { x: 0., y: 0. }; 100]); 100]),
            resolution: Resolution {
                width: 100,
                height: 100,
            },
            cell_size: 10.,
        }
    }
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_precision_loss,
        clippy::cast_sign_loss
    )]
    fn process_xy(&self, x: f32, y: f32) -> Vec2 {
        let mut x = ((x + X_EXTENT) / self.cell_size) as usize;
        let mut y = ((y + Y_EXTENT) / self.cell_size) as usize;
        if x >= self.resolution.width {
            x = self.resolution.width - 1;
        }
        if y >= self.resolution.height {
            y = self.resolution.height - 1;
        }
        self.grid[y][x]
    }
    fn collect_peak(&mut self, peak_x: f32, peak_y: f32, peak_amplitude: f32) {
        self.peaks.push(Vec2::new(peak_x, peak_y));
        self.peak_amplitudes.push(-peak_amplitude);
    }
    fn clear_peaks(&mut self) {
        self.peaks.clear();
        self.peak_amplitudes.clear();
        self.grid = Grid(vec![SubGrid(vec![Vec2 { x: 0., y: 0. }; 100]); 100]);
    }
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_precision_loss,
        clippy::cast_sign_loss
    )]
    fn compute_influence(&mut self) {
        for (i, peak) in self.peaks.iter().enumerate() {
            for x in 0..self.resolution.width {
                let x_coord = (x as f32).mul_add(self.cell_size, -X_EXTENT);
                for y in 0..self.resolution.height {
                    let y_coord = (y as f32).mul_add(self.cell_size, -Y_EXTENT);
                    let distance = vec2(x_coord + 5., y_coord + 5.).distance(*peak);
                    let amplitude = self.peak_amplitudes[i] * 0.01;
                    let influence = amplitude * (-distance.powi(2) / (2. * self.width.powi(2))).exp();
                    let x_influence = (x_coord - peak.x) * influence;
                    let y_influence = (y_coord - peak.y) * influence;
                    self.grid[y][x] += vec2(x_influence, y_influence);
                }
            }
        }
    }
}
#[allow(clippy::needless_pass_by_value)]
fn move_shapes(
    mut query: Query<(&mut Transform, &mut Velocity, &Mesh2dHandle), With<Velocity>>,
    mut system: Query<&mut System>,
) {
    for mut entity in &mut query {
        let delta = system
            .single_mut()
            .process_xy(entity.0.translation.x, entity.0.translation.y);
        let mass = entity.0.scale.x;
        entity.1.speed_x += delta.x / mass;
        entity.1.speed_y += delta.y / mass;
        // Slightly draw each shape back to the center of the screen
        entity.1.speed_x -= entity.0.translation.x * 0.003;
        entity.1.speed_y -= entity.0.translation.y * 0.003;
        // Apply a slight damping to the speed
        entity.1.speed_x *= 0.9999;
        entity.1.speed_y *= 0.9999;
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
    for x in -10..10 {
        for y in -10..10 {
            let translation = Vec3::new(x as f32 * 50., y as f32 * 50., -10.);
            let grid_velocity = system.single().process_xy(translation.x, translation.y);
            let radius = grid_velocity.length().mul_add(10000., 5.);
            let dot_radius = std::ops::Mul::mul(radius, 0.0005).clamp(0.5, 20.);
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
