import pygame
import numpy as np

# Конфигурация
WINDOW_WIDTH = 1400
WINDOW_HEIGHT = 720
BACKGROUND_COLOR = (30, 30, 30)
PARTICLE_COLOR = (100, 200, 255)
PARTICLE_RADIUS = 2
SIM_SCALE = 0.7  # Масштаб подберите под ваш случай
SIM_OFFSET_X = WINDOW_WIDTH // 2  # Центрирование по X
SIM_OFFSET_Y = 100  # Отступ сверху


def load_frames(filename):
    """Загрузка всех кадров из файла"""
    frames = []
    with open(filename, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            data = np.fromstring(line, sep=' ', dtype=np.float32)
            if len(data) % 2 != 0:
                continue  # Пропустить некорректные строки
            frames.append(data.reshape(-1, 2))
    return frames


def main():
    pygame.init()
    screen = pygame.display.set_mode((WINDOW_WIDTH, WINDOW_HEIGHT))
    pygame.display.set_caption("Fluid Simulation Visualizer")
    clock = pygame.time.Clock()

    frames = load_frames("OutputData/res.txt")
    if not frames:
        print("No frames loaded!")
        return

    current_frame = 0
    paused = False
    running = True

    while running:
        # Обработка событий
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    paused = not paused
                elif event.key == pygame.K_LEFT:
                    current_frame = max(0, current_frame - 1)
                elif event.key == pygame.K_RIGHT:
                    current_frame = min(len(frames) - 1, current_frame + 1)
                elif event.key == pygame.K_ESCAPE:
                    running = False

        # Отрисовка
        screen.fill(BACKGROUND_COLOR)

        if not paused and current_frame < len(frames):
            particles = frames[current_frame]

            # Рисуем все частицы
            for x, y in particles:
                screen_x = int(x * SIM_SCALE)
                screen_y = int(WINDOW_HEIGHT - (y * SIM_SCALE + SIM_OFFSET_Y))

                pygame.draw.circle(
                    screen,
                    PARTICLE_COLOR,
                    (screen_x, screen_y),
                    PARTICLE_RADIUS
                )

            current_frame += 1

        pygame.display.flip()
        clock.tick(60)  # Максимум 60 FPS

    pygame.quit()


if __name__ == "__main__":
    main()