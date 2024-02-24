use bevy::{
    math::vec2,
    prelude::*,
    render::{
        settings::{Backends, RenderCreation, WgpuSettings},
        RenderPlugin,
    },
    sprite::{MaterialMesh2dBundle, Mesh2dHandle},
};
const X_EXTENT: f32 = 500.;
const Y_EXTENT: f32 = 500.;
const DIMENSION: usize = 100;
const CELL_SIZE: f32 = 10.;
#[derive(Component)]
struct Velocity {
    speed_x: f32,
    speed_y: f32,
}
#[derive(Component)]
struct GridVelocity;
#[derive(Clone)]
struct Grid {
    data: Vec<Vec<f32>>,
}
impl Grid {
    fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![vec![0.; width]; height],
        }
    }
    fn diffuse(&mut self, diffusion_rate: f32) {
        // Assume the grid is square
        let width = self.data.len();
        if width == 0 {
            return;
        }
        let divisor = 4.0_f32.mul_add(diffusion_rate, 1.0);
        for y in 1..width - 1 {
            for x in 1..width - 1 {
                self.data[y][x] = (self.data[y][x - 1]
                    + self.data[y][x + 1]
                    + self.data[y - 1][x]
                    + self.data[y + 1][x])
                    .mul_add(diffusion_rate, self.data[y][x])
                    / divisor;
            }
        }
        let edge_rate = 0.75;
        for x in 0..width {
            self.data[x][0] *= edge_rate;
            self.data[x][width - 1] *= edge_rate;
            self.data[0][x] *= edge_rate;
            self.data[width - 1][x] *= edge_rate;
        }
    }
}
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
            speed_x: 10.,
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
    let spawn_dots = true;
    if spawn_dots {
        for _ in 0..10 {
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
    }
    commands.spawn(System::new());
}
#[derive(Component)]
struct System {
    peaks: Vec<Matter>,
    grid: Grid,
    resolution: Resolution, // Defines the resolution of the grid (width, height)
    cell_size: f32,         // The size of each cell in the grid
}
struct Matter {
    grid_x: usize,
    grid_y: usize,
    mass: f32,
}
struct Resolution {
    width: usize,
    height: usize,
}
impl System {
    fn new() -> Self {
        Self {
            peaks: Vec::new(),
            grid: Grid::new(DIMENSION, DIMENSION),
            resolution: Resolution {
                width: DIMENSION,
                height: DIMENSION,
            },
            cell_size: CELL_SIZE,
        }
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
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
        self.peaks.push(Matter {
            grid_x: self.position_to_cell(peak_x),
            grid_y: self.position_to_cell(peak_y),
            mass: peak_amplitude,
        });
    }
    fn compute_influence(&mut self) {
        for peak in &self.peaks {
            self.grid.data[peak.grid_y][peak.grid_x] += peak.mass;
        }
    }
    fn diffuse(&mut self, steps: usize, diffusion_rate: f32) {
        for _ in 0..steps {
            self.grid.diffuse(diffusion_rate);
        }
    }
}
#[allow(clippy::needless_pass_by_value)]
fn move_shapes(
    time: Res<Time>,
    mut query: Query<(&mut Transform, &mut Velocity, &Mesh2dHandle), With<Velocity>>,
    system: Query<&mut System>,
) {
    let time_delta = time.delta_seconds();
    for mut entity in &mut query {
        let mut delta = vec2(0., 0.);
        let cell_x = system.single().position_to_cell(entity.0.translation.x);
        let cell_y = system.single().position_to_cell(entity.0.translation.y);
        if cell_x < 1
            || cell_x >= system.single().resolution.width - 1
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
        let left_influence = system.single().grid.data[cell_y][cell_x - 1];
        let right_influence = system.single().grid.data[cell_y][cell_x + 1];
        let up_influence = system.single().grid.data[cell_y - 1][cell_x];
        let down_influence = system.single().grid.data[cell_y + 1][cell_x];
        delta.x = right_influence - left_influence;
        delta.y = down_influence - up_influence;
        let mass = entity.0.scale.x;
        if mass > 0. {
            entity.1.speed_x += delta.x / mass;
            entity.1.speed_y += delta.y / mass;
        }
        let lensing_strength = 0.1;
        entity.1.speed_x += delta.x * lensing_strength;
        entity.1.speed_y += delta.y * lensing_strength;

        entity.0.translation.x += entity.1.speed_x * time_delta;
        entity.0.translation.y += entity.1.speed_y * time_delta;
    }
}
#[allow(clippy::needless_pass_by_value)]
fn add_peaks(
    mut system: Query<&mut System>,
    query: Query<(&Transform, &Mesh2dHandle), With<Velocity>>,
) {
    system.single_mut().peaks.clear();
    for entity in &mut query.iter() {
        let radius = entity.0.scale.x;
        system
            .single_mut()
            .collect_peak(entity.0.translation.x, entity.0.translation.y, radius);
    }
    system.single_mut().compute_influence();
}
#[allow(
    clippy::needless_pass_by_value,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation
)]
fn diffuse(time: Res<Time>, mut system: Query<&mut System>) {
    let time_delta = time.delta_seconds();
    let steps = (time_delta * 60.).round() as usize;
    system.single_mut().diffuse(steps, 0.8);
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
    for x in 0..system.single().resolution.width {
        if x % 4 != 0 {
            continue;
        }
        for y in 0..system.single().resolution.height {
            if y % 4 != 0 {
                continue;
            }
            let translation = Vec3::new(
                (x as f32 - system.single().resolution.width as f32 / 2.)
                    .mul_add(system.single().cell_size, system.single().cell_size / 2.),
                (y as f32 - system.single().resolution.height as f32 / 2.)
                    .mul_add(system.single().cell_size, system.single().cell_size / 2.),
                -10.,
            );
            let grid_velocity = system.single().grid.data[y][x];
            let radius = grid_velocity.mul_add(10., 5.);
            let dot_radius = radius.abs().clamp(0.2, 20.);
            let dot = Mesh2dHandle(meshes.add(Circle { radius: dot_radius }));
            let hue = (grid_velocity + 1.) / 2.;
            let hue = hue * 360.;
            let saturation = // out of 1.0
                1. - (grid_velocity / 100.).min(0.01);
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
        .add_systems(Update, diffuse.after(add_peaks).before(move_shapes))
        .add_systems(Update, change_dot_color.after(add_peaks))
        .run();
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_position_to_cell() {
        let system = crate::System::new();
        let dimension_half = DIMENSION / 2;
        assert_eq!(system.position_to_cell(-500.), 0);
        assert_eq!(system.position_to_cell(500.), DIMENSION);
        assert_eq!(system.position_to_cell(0.), dimension_half);
    }
    #[test]
    fn test_collect_peak() {
        let mut system = crate::System::new();
        system.collect_peak(0., 0., 1.);
        assert_eq!(system.peaks.len(), 1);
        assert_eq!(system.peaks[0].grid_x, 50);
        assert_eq!(system.peaks[0].grid_y, 50);
        let error_margin = f32::EPSILON;
        assert!((system.peaks[0].mass - 1.).abs() < error_margin);
    }
    #[test]
    fn test_compute_influence() {
        let mut system = crate::System::new();
        system.collect_peak(0., 0., 1.);
        system.compute_influence();
        let error_margin = f32::EPSILON;
        assert!((system.grid.data[50][50] - 1.) < error_margin);
    }
}
