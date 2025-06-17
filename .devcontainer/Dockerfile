# C++ Dev Container Dockerfile
FROM mcr.microsoft.com/devcontainers/cpp:ubuntu-24.04

# Install additional packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gdb \
        clang-format \
        zsh \
        git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set up a non-root user (optional, VS Code dev containers best practice)
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=1000
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && apt-get update && apt-get install -y sudo \
    && echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

USER $USERNAME

# Install zsh plugins and set up p10k for the non-root user
RUN git clone --depth=1 https://github.com/romkatv/powerlevel10k.git /home/$USERNAME/powerlevel10k \
    && git clone https://github.com/zsh-users/zsh-autosuggestions /home/$USERNAME/.zsh-autosuggestions \
    && echo "source /home/$USERNAME/powerlevel10k/powerlevel10k.zsh-theme" >> /home/$USERNAME/.zshrc \
    && echo "source /home/$USERNAME/.zsh-autosuggestions/zsh-autosuggestions.zsh" >> /home/$USERNAME/.zshrc \
    && chsh -s $(which zsh) $USERNAME

WORKDIR /workspace
