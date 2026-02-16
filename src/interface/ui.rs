use ratatui::{
    prelude::*,
    widgets::{
        Block, Borders, BorderType, Paragraph, Tabs, Gauge, 
        Sparkline, Table, Row, Cell, Wrap, ListItem, List,
        canvas::{Canvas, Circle, Line as CanvasLine},
    },
    style::{Color, Style, Modifier},
    text::{Line, Span},
};
use crate::interface::state::{AppState, AppMode, WorkerStatus};
use crate::core::domain::Cluster;

// --- Color Palette ---
const COL_BG: Color = Color::Reset;
const COL_FG: Color = Color::White;
const COL_HIGHLIGHT: Color = Color::Yellow;
const COL_ACCENT: Color = Color::Cyan;
const COL_ATOM_A: Color = Color::LightRed;   // Oxygen / Anion
const COL_ATOM_B: Color = Color::LightBlue;  // Magnesium / Cation
const COL_BOND: Color = Color::DarkGray;
const COL_SUCCESS: Color = Color::Green;
const COL_FAIL: Color = Color::Red;
const COL_HEADER: Color = Color::Magenta;
const COL_DIVERSITY: Color = Color::LightGreen;

pub fn draw(f: &mut Frame, app: &mut AppState) {
    if f.area().width < 40 || f.area().height < 10 {
        let p = Paragraph::new("Terminal too small.")
            .alignment(Alignment::Center)
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(p, f.area());
        return;
    }

    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3), 
            Constraint::Min(0),    
            Constraint::Length(1), 
        ])
        .split(f.area());

    draw_header(f, app, chunks[0]);
    
    match app.mode {
        AppMode::Dashboard => draw_dashboard(f, app, chunks[1]),
        AppMode::HallOfFame => draw_hall_of_fame(f, app, chunks[1]),
        AppMode::Analysis => draw_analysis(f, app, chunks[1]),
        AppMode::StructureViewer => draw_fullscreen_viewer(f, app, chunks[1]),
        _ => {}
    }

    draw_footer(f, app, chunks[2]);
}

fn draw_header(f: &mut Frame, app: &AppState, area: Rect) {
    let titles = vec![" 1:Dash ", " 2:Analysis ", " 3:HoF ", " 4:Viewer "];
    let idx = match app.mode {
        AppMode::Dashboard => 0,
        AppMode::Analysis => 1,
        AppMode::HallOfFame => 2,
        AppMode::StructureViewer => 3,
        _ => 0,
    };

    let tabs = Tabs::new(titles)
        .block(Block::default().borders(Borders::BOTTOM))
        .select(idx)
        .highlight_style(Style::default().fg(COL_HIGHLIGHT).add_modifier(Modifier::BOLD));
    
    f.render_widget(tabs, area);
}

fn draw_footer(f: &mut Frame, app: &AppState, area: Rect) {
    let status_str = match app.worker_status {
        WorkerStatus::Running => "RUNNING",
        WorkerStatus::Paused => "PAUSED",
        WorkerStatus::Idle => "IDLE",
        WorkerStatus::Starting => "STARTING",
        WorkerStatus::Finished => "DONE",
        WorkerStatus::Error => "ERROR",
    };
    
    let color = match app.worker_status {
        WorkerStatus::Running => COL_SUCCESS,
        WorkerStatus::Error => COL_FAIL,
        WorkerStatus::Paused => COL_HIGHLIGHT,
        _ => COL_FG,
    };

    let best_val = app.current_best.as_ref()
        .and_then(|c| c.energy)
        .unwrap_or(0.0);

    let text = Line::from(vec![
        Span::styled(format!(" STATUS: {:<8}", status_str), Style::default().fg(color).add_modifier(Modifier::BOLD)),
        Span::raw(" | "),
        Span::raw(format!("Ops/s: {:<6.1}", app.ops_per_second)),
        Span::raw(" | "),
        Span::styled(format!("Best: {:.4} eV", best_val), Style::default().fg(COL_ACCENT)),
        Span::raw(" | [Q]uit [Space]Pause [R]eset-View"),
    ]);

    let p = Paragraph::new(text)
        .style(Style::default().bg(Color::DarkGray).fg(Color::White));
    f.render_widget(p, area);
}

fn draw_dashboard(f: &mut Frame, app: &AppState, area: Rect) {
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(60), Constraint::Percentage(40)])
        .split(area);

    let left_rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(65), Constraint::Percentage(35)])
        .split(cols[0]);

    if let Some(cluster) = &app.active_cluster {
        draw_cluster_3d(f, app, left_rows[0], cluster, " Live Structure ");
    } else {
        f.render_widget(Block::default().title(" Waiting for Data... ").borders(Borders::ALL), left_rows[0]);
    }

    draw_population_charts(f, app, left_rows[1]);

    let right_rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(50), Constraint::Percentage(20), Constraint::Percentage(30)])
        .split(cols[1]);

    draw_logs(f, app, right_rows[0]);
    draw_gauges(f, app, right_rows[1]);
    draw_stats(f, app, right_rows[2]);
}

fn draw_cluster_3d(f: &mut Frame, app: &AppState, area: Rect, cluster: &Cluster, title: &str) {
    let block = Block::default()
        .title(title)
        .borders(Borders::ALL)
        .border_type(BorderType::Rounded);
    
    let inner_area = block.inner(area);
    f.render_widget(block, area);

    if inner_area.width < 1 || inner_area.height < 1 { return; }
    if cluster.atoms.is_empty() { return; }

    let mut render_atoms: Vec<(f64, f64, f64, Color, f64)> = cluster.atoms.iter().map(|a| {
        use nalgebra::{Rotation3, Vector3};
        let rot_y = Rotation3::from_axis_angle(&Vector3::y_axis(), app.viewport.azimuth);
        let rot_x = Rotation3::from_axis_angle(&Vector3::x_axis(), app.viewport.elevation);
        let p_rot = rot_x * rot_y * a.position;
        
        let color = if a.element_id == 0 { COL_ATOM_A } else { COL_ATOM_B };
        let size = if a.element_id == 0 { 0.6 } else { 0.4 }; 
        (p_rot.x, p_rot.y, p_rot.z, color, size)
    }).collect();

    if render_atoms.iter().any(|(x, y, z, _, _)| x.is_nan() || y.is_nan() || z.is_nan()) {
        f.render_widget(Paragraph::new("Error: NaN Coordinates").style(Style::default().fg(COL_FAIL)), inner_area);
        return;
    }

    let max_coord = render_atoms.iter()
        .flat_map(|(x, y, _, _, _)| vec![x.abs(), y.abs()])
        .fold(0.0, f64::max)
        .max(1.0); 
    
    let bound = max_coord * 1.2;

    render_atoms.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

    let canvas = Canvas::default()
        .background_color(COL_BG)
        .x_bounds([-bound, bound])
        .y_bounds([-bound, bound])
        .paint(|ctx| {
            for i in 0..render_atoms.len() {
                for j in (i+1)..render_atoms.len() {
                    let a = &render_atoms[i];
                    let b = &render_atoms[j];
                    let dx = a.0 - b.0;
                    let dy = a.1 - b.1;
                    let dz = a.2 - b.2;
                    let d2 = dx*dx + dy*dy + dz*dz;
                    
                    if d2 < 5.0 { 
                        let avg_z = (a.2 + b.2) / 2.0;
                        let brightness = if avg_z < 0.0 { Color::DarkGray } else { COL_BOND };
                        ctx.draw(&CanvasLine {
                            x1: a.0, y1: a.1,
                            x2: b.0, y2: b.1,
                            color: brightness,
                        });
                    }
                }
            }
            for (x, y, z, col, size) in &render_atoms {
                let perspective = (1.0 + z * 0.05).clamp(0.5, 1.5);
                let view_zoom = app.viewport.zoom;
                ctx.draw(&Circle {
                    x: *x * view_zoom,
                    y: *y * view_zoom,
                    radius: *size * view_zoom * perspective,
                    color: *col,
                });
            }
        });
    
    f.render_widget(canvas, inner_area);

    let rot_status = if app.viewport.auto_rotate { "Auto-Rot: ON" } else { "Auto-Rot: OFF" };
    let overlay = Paragraph::new(rot_status).style(Style::default().fg(Color::DarkGray).add_modifier(Modifier::ITALIC));
    let overlay_area = Rect { x: inner_area.x + inner_area.width.saturating_sub(14), y: inner_area.y, width: 14, height: 1 };
    f.render_widget(overlay, overlay_area);
}

fn draw_population_charts(f: &mut Frame, app: &AppState, area: Rect) {
    let block = Block::default().title(" Population Analysis ").borders(Borders::ALL);
    let inner = block.inner(area);
    f.render_widget(block, area);

    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(50), Constraint::Percentage(50)])
        .split(inner);

    if !app.telemetry.best_energy_history.is_empty() {
        // FIXED: Use correct field names
        let min = app.telemetry.global_min_energy;
        let max = app.telemetry.global_max_energy;
        let range = (max - min).max(0.1);
        
        let width = inner.width as usize;
        let data: Vec<u64> = app.telemetry.best_energy_history.iter()
            .rev()
            .take(width)
            .map(|(_, e)| {
                let norm = (e - min) / range;
                (norm * 10.0) as u64
            })
            .collect();
        
        let data_rev: Vec<u64> = data.into_iter().rev().collect();

        let spark = Sparkline::default()
            .block(Block::default().title("Energy Convergence").borders(Borders::NONE))
            .style(Style::default().fg(COL_ACCENT))
            .data(&data_rev);
        f.render_widget(spark, chunks[0]);
    }

    if !app.telemetry.diversity_history.is_empty() {
        let width = inner.width as usize;
        let data: Vec<u64> = app.telemetry.diversity_history.iter()
            .rev()
            .take(width)
            .map(|(_, d)| {
                (d / 10.0) as u64
            })
            .collect();
        let data_rev: Vec<u64> = data.into_iter().rev().collect();

        let spark_div = Sparkline::default()
            .block(Block::default().title("Genetic Diversity").borders(Borders::NONE))
            .style(Style::default().fg(COL_DIVERSITY))
            .data(&data_rev);
        f.render_widget(spark_div, chunks[1]);
    }
}

fn draw_logs(f: &mut Frame, app: &AppState, area: Rect) {
    let block = Block::default().title(" System Log ").borders(Borders::ALL);
    let inner = block.inner(area);
    f.render_widget(block, area);

    let items: Vec<ListItem> = app.logs.iter().rev().map(|line| {
        let style = if line.to_lowercase().contains("error") || line.contains("failed") {
            Style::default().fg(COL_FAIL)
        } else if line.contains(">>>") || line.contains("Extinction") {
            Style::default().fg(COL_SUCCESS)
        } else {
            Style::default().fg(Color::Gray)
        };

        ListItem::new(Line::from(vec![
            Span::styled(">", Style::default().fg(Color::DarkGray)),
            Span::raw(" "),
            Span::raw(line),
        ])).style(style)
    }).collect();

    let list = List::new(items);
    f.render_widget(list, inner);
}

fn draw_gauges(f: &mut Frame, app: &AppState, area: Rect) {
    let block = Block::default().title(" Health ").borders(Borders::ALL);
    let inner = block.inner(area);
    f.render_widget(block, area);

    let layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(1); 2])
        .split(inner);

    let div = app.telemetry.diversity_history.back().map(|&(_, d)| d).unwrap_or(100.0);
    let ratio = div / 100.0;
    
    let label = format!("Diversity: {:.0}%", div);
    let color = if ratio > 0.4 { COL_DIVERSITY } else { COL_FAIL };
    
    let g_div = Gauge::default()
        .gauge_style(Style::default().fg(color).bg(Color::DarkGray))
        .ratio(ratio.clamp(0.0, 1.0))
        .label(label);
    
    f.render_widget(g_div, layout[0]);
}

fn draw_stats(f: &mut Frame, app: &AppState, area: Rect) {
    let block = Block::default().title(" Statistics ").borders(Borders::ALL);
    let inner = block.inner(area);
    f.render_widget(block, area);

    // FIXED: Use correct global_min_energy field
    let text = vec![
        Line::from(vec![
            Span::styled("Algorithm: ", Style::default().fg(Color::Gray)),
            Span::styled(format!("{:?}", app.params.algorithm), Style::default().fg(COL_HIGHLIGHT))
        ]),
        Line::from(vec![
            Span::styled("Iteration: ", Style::default().fg(Color::Gray)),
            Span::styled(format!("{}", app.total_iterations), Style::default().fg(COL_HIGHLIGHT))
        ]),
        Line::from(vec![
            Span::styled("Elite Count: ", Style::default().fg(Color::Gray)),
            Span::styled(format!("{}", app.hall_of_fame.len()), Style::default().fg(COL_HIGHLIGHT))
        ]),
        Line::from(vec![
            Span::styled("Global Min: ", Style::default().fg(Color::Gray)),
            Span::styled(format!("{:.5} eV", app.telemetry.global_min_energy), Style::default().fg(COL_SUCCESS).add_modifier(Modifier::BOLD))
        ]),
    ];

    f.render_widget(Paragraph::new(text).wrap(Wrap { trim: true }), inner);
}

fn draw_hall_of_fame(f: &mut Frame, app: &mut AppState, area: Rect) {
    let header_cells = ["Rank", "ID", "Energy (eV)", "Origin", "Hash"]
        .iter()
        .map(|h| Cell::from(*h).style(Style::default().fg(COL_HEADER)));
    let header = Row::new(header_cells).height(1).bottom_margin(1);

    if app.hall_of_fame.is_empty() {
        f.render_widget(
            Paragraph::new("Waiting for successful evaluations...")
                .block(Block::default().borders(Borders::ALL).title(" Hall of Fame "))
                .alignment(Alignment::Center), 
            area
        );
        return;
    }

    let rows = app.hall_of_fame.iter().enumerate().map(|(i, c)| {
        let id_short = c.id.to_string();
        let id_disp = if id_short.len() > 8 { &id_short[..8] } else { &id_short };
        let hash_disp = c.hash_key.as_deref().unwrap_or("-");
        let hash_short = if hash_disp.len() > 12 { &hash_disp[..12] } else { hash_disp };

        let cells = vec![
            Cell::from(format!("#{}", i + 1)),
            Cell::from(id_disp.to_string()),
            Cell::from(format!("{:.5}", c.energy.unwrap_or(0.0))),
            Cell::from(c.origin.clone()),
            Cell::from(hash_short.to_string()),
        ];
        Row::new(cells).height(1)
    });

    let t = Table::new(rows, &[
        Constraint::Length(6),
        Constraint::Length(10),
        Constraint::Length(15),
        Constraint::Length(15),
        Constraint::Min(10),
    ])
    .header(header)
    .block(Block::default().borders(Borders::ALL).title(format!(" Hall of Fame ({}) ", app.hall_of_fame.len())))
    .row_highlight_style(Style::default().add_modifier(Modifier::REVERSED));

    f.render_stateful_widget(t, area, &mut app.hof_state);
}

fn draw_analysis(f: &mut Frame, app: &AppState, area: Rect) {
    draw_config(f, app, area);
}

fn draw_config(f: &mut Frame, app: &AppState, area: Rect) {
    let block = Block::default().borders(Borders::ALL).title(" Simulation Parameters ");
    let inner = block.inner(area);
    f.render_widget(block, area);

    let p = &app.params;
    
    let kv = |k: &str, v: String| -> ListItem {
        ListItem::new(Line::from(vec![
            Span::styled(format!("{:<15}", k), Style::default().fg(COL_ACCENT)), 
            Span::raw(v)
        ]))
    };

    let items = vec![
        kv("Atom Count:", p.atom_count.to_string()),
        kv("Box Size:", format!("{:.1} Å", p.box_size)),
        kv("Threads:", p.threads.to_string()),
        ListItem::new(Line::from(" ")),
        kv("Algorithm:", format!("{:?}", p.algorithm)),
        kv("Pop Size:", p.population_size.to_string()),
        kv("Mutation Rate:", format!("{:.2}", p.mutation_rate)),
        kv("Crossover:", format!("{:.2}", p.crossover_rate)),
        ListItem::new(Line::from(" ")),
        kv("Temperature:", format!("{:.1} K", p.temperature)),
        kv("Step Size:", format!("{:.2} Å", p.step_size)),
    ];

    let list = List::new(items).block(Block::default().borders(Borders::NONE));
    f.render_widget(list, inner);
}

fn draw_fullscreen_viewer(f: &mut Frame, app: &AppState, area: Rect) {
    if let Some(cluster) = &app.active_cluster {
        draw_cluster_3d(f, app, area, cluster, " Structure Viewer (Fullscreen) ");
    } else {
        let p = Paragraph::new("No structure selected.\nGo to Hall of Fame and select a structure.")
            .alignment(Alignment::Center)
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(p, area);
    }
}