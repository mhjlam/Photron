#include "window.hpp"

bool GLWindow::create(const std::string& title, int width, int height) {
	// Register window class
	WNDCLASSEXA wc = {};
	wc.cbSize = sizeof(WNDCLASSEXA);
	wc.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
	wc.lpfnWndProc = window_proc;
	wc.hInstance = GetModuleHandle(nullptr);
	wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
	wc.lpszClassName = "PhotonGLWindow";

	if (!RegisterClassExA(&wc)) {
		std::cerr << "Failed to register window class" << std::endl;
		return false;
	}

	// Create window
	hwnd_ = CreateWindowExA(0, "PhotonGLWindow", title.c_str(), WS_OVERLAPPEDWINDOW | WS_VISIBLE, CW_USEDEFAULT,
							CW_USEDEFAULT, width, height, nullptr, nullptr, GetModuleHandle(nullptr), this);

	if (!hwnd_) {
		std::cerr << "Failed to create window" << std::endl;
		return false;
	}

	// Get device context
	hdc_ = GetDC(hwnd_);
	if (!hdc_) {
		std::cerr << "Failed to get device context" << std::endl;
		return false;
	}

	// Choose pixel format
	PIXELFORMATDESCRIPTOR pfd = {};
	pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR);
	pfd.nVersion = 1;
	pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
	pfd.iPixelType = PFD_TYPE_RGBA;
	pfd.cColorBits = 32;
	pfd.cDepthBits = 24;
	pfd.cStencilBits = 8;

	int pixel_format = ChoosePixelFormat(hdc_, &pfd);
	if (!pixel_format) {
		std::cerr << "Failed to choose pixel format" << std::endl;
		return false;
	}

	if (!SetPixelFormat(hdc_, pixel_format, &pfd)) {
		std::cerr << "Failed to set pixel format" << std::endl;
		return false;
	}

	// Create OpenGL context
	hglrc_ = wglCreateContext(hdc_);
	if (!hglrc_) {
		std::cerr << "Failed to create OpenGL context" << std::endl;
		return false;
	}

	if (!wglMakeCurrent(hdc_, hglrc_)) {
		std::cerr << "Failed to make OpenGL context current" << std::endl;
		return false;
	}

	// Print OpenGL info
	std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;
	std::cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << std::endl;
	std::cout << "OpenGL Vendor: " << glGetString(GL_VENDOR) << std::endl;

	return true;
}

void GLWindow::swap_buffers() {
	if (hdc_) {
		SwapBuffers(hdc_);
	}
}

bool GLWindow::should_close() {
	return should_close_;
}

void GLWindow::poll_events() {
	// Reset mouse deltas at the start of each frame
	mouse_state_.reset_deltas();

	MSG msg;
	while (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)) {
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
}

void GLWindow::cleanup() {
	if (hglrc_) {
		wglMakeCurrent(nullptr, nullptr);
		wglDeleteContext(hglrc_);
		hglrc_ = nullptr;
	}

	if (hdc_) {
		ReleaseDC(hwnd_, hdc_);
		hdc_ = nullptr;
	}

	if (hwnd_) {
		DestroyWindow(hwnd_);
		hwnd_ = nullptr;
	}
}

LRESULT CALLBACK GLWindow::window_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam) {
	GLWindow* window = nullptr;

	if (msg == WM_NCCREATE) {
		CREATESTRUCT* create = reinterpret_cast<CREATESTRUCT*>(lparam);
		window = static_cast<GLWindow*>(create->lpCreateParams);
		SetWindowLongPtr(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(window));
	}
	else {
		window = reinterpret_cast<GLWindow*>(GetWindowLongPtr(hwnd, GWLP_USERDATA));
	}

	if (window) {
		switch (msg) {
			case WM_CLOSE: {
				window->should_close_ = true;
				return 0;
			}
			case WM_KEYDOWN: {
				if (wparam == VK_ESCAPE) {
					window->should_close_ = true;
				}
				return 0;
			}

			// Mouse button events
			case WM_LBUTTONDOWN: {
				window->mouse_state_.left_button_down = true;
				SetCapture(hwnd);
				return 0;
			}
			case WM_LBUTTONUP: {
				window->mouse_state_.left_button_down = false;
				ReleaseCapture();
				return 0;
			}
			case WM_RBUTTONDOWN: {
				window->mouse_state_.right_button_down = true;
				SetCapture(hwnd);
				return 0;
			}
			case WM_RBUTTONUP: {
				window->mouse_state_.right_button_down = false;
				ReleaseCapture();
				return 0;
			}
			case WM_MBUTTONDOWN: {
				window->mouse_state_.middle_button_down = true;
				SetCapture(hwnd);
				return 0;
			}
			case WM_MBUTTONUP: {
				window->mouse_state_.middle_button_down = false;
				ReleaseCapture();
				return 0;
			}

			// Mouse movement
			case WM_MOUSEMOVE: {
				int x = LOWORD(lparam);
				int y = HIWORD(lparam);
				window->mouse_state_.update(x, y);
				return 0;
			}

			// Mouse wheel
			case WM_MOUSEWHEEL: {
				short wheel_delta = GET_WHEEL_DELTA_WPARAM(wparam);
				window->mouse_state_.scroll_delta = static_cast<float>(wheel_delta) / WHEEL_DELTA;
				return 0;
			}
		}
	}

	return DefWindowProc(hwnd, msg, wparam, lparam);
}
