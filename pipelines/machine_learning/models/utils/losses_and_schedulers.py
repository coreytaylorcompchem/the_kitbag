import torch
import torch.nn as nn
import torch.optim as optim

class LabelSmoothingLoss(nn.Module):
    def __init__(self, classes, smoothing=0.1, ignore_index=-100):
        super().__init__()
        self.confidence = 1.0 - smoothing
        self.smoothing = smoothing
        self.cls = classes
        self.ignore_index = ignore_index

    def forward(self, pred, target):
        pred = pred.log_softmax(dim=-1)
        with torch.no_grad():
            true_dist = torch.zeros_like(pred)
            true_dist.fill_(self.smoothing / (self.cls - 1))
            mask = target == self.ignore_index
            target = target.masked_fill(mask, 0)
            true_dist.scatter_(1, target.unsqueeze(1), self.confidence)
            true_dist.masked_fill_(mask.unsqueeze(1), 0)
        return torch.mean(torch.sum(-true_dist * pred, dim=-1))


class WarmupInverseSqrtScheduler(optim.lr_scheduler._LRScheduler):
    def __init__(self, optimizer, warmup_steps, last_epoch=-1):
        self.warmup_steps = warmup_steps
        super().__init__(optimizer, last_epoch)

    def get_lr(self):
        step = max(self.last_epoch, 1)
        scale = (self.warmup_steps ** 0.5) * min(step ** (-0.5), step * (self.warmup_steps ** (-1.5)))
        return [base_lr * scale for base_lr in self.base_lrs]
